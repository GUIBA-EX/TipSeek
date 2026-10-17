use ahash::AHashMap;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::env;
use std::ffi::{c_char, c_int, c_uint, c_void, CString};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

const CACHE_MAGIC: &[u8; 4] = b"GM2K";
const CACHE_VERSION: u16 = 4;
const FLAG_CANONICAL_KMERS: u16 = 1;
const MEBIBYTE_READS: u64 = 1_048_576;
const TOTAL_BUFFER_BUDGET: usize = 64 * 1024 * 1024;
const BUFFER_LOW_WATERMARK: usize = TOTAL_BUFFER_BUDGET / 2;
const MIN_FILE_BUDGET: usize = 128 * 1024;
// 输入端也留一根粗管子:少喊几次 gzread/read,解压和系统调用的开销摊得更薄。
const READ_BUFFER_SIZE: usize = 1024 * 1024;

type AppResult<T> = Result<T, String>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FileKind {
    Fasta,
    Fastq,
}

#[derive(Clone, Debug)]
struct Args {
    reference: PathBuf,
    q1: Vec<PathBuf>,
    q2: Vec<PathBuf>,
    kmer: usize,
    step: usize,
    max_read_blocks: u64,
    output: PathBuf,
    out_subdir: String,
    dictionary: Option<PathBuf>,
    get_reverse: bool,
    use_composition_pattern: bool,
    mode: u8,
    fallback_kmers: Vec<usize>,
    minimum_hit_span: usize,
    link_rad_arms: bool,
    link_rad_max_fragments: u64,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            reference: PathBuf::new(),
            q1: Vec::new(),
            q2: Vec::new(),
            kmer: 31,
            step: 3,
            max_read_blocks: 0,
            output: PathBuf::new(),
            out_subdir: "filtered".to_string(),
            dictionary: None,
            get_reverse: false,
            use_composition_pattern: false,
            mode: 0,
            fallback_kmers: Vec::new(),
            minimum_hit_span: 0,
            link_rad_arms: false,
            link_rad_max_fragments: 256,
        }
    }
}

// 把能用的开关一次说明白，用户别靠猜参数整活儿。
fn print_help() {
    println!(
        "MainFilterNew (Rust implementation)\n\
         Usage: MainFilterNew -r REF -o OUT [options]\n\n\
         -r PATH       Reference FASTA file or directory\n\
         -q1 FILE...   Read 1 / single-end FASTA, FASTQ, or gzipped FASTQ\n\
         -q2 FILE...   Read 2 files\n\
         -kf INT       K-mer length (default: 31)\n\
         -s INT        Read scanning step (default: 3)\n\
         -m_reads INT  Maximum 2^20 read records per input file; 0 disables\n\
         -o PATH       Output directory\n\
         -subdir NAME  Filtered-read subdirectory (default: filtered)\n\
         -lkd PATH     Load or write a reusable k-mer dictionary\n\
         -gr           Index reverse-complement reference k-mers\n\
         -lb           Accept the legacy composition-pattern option\n\
         -m INT        Mode 0..5 (default: 0)\n\
         --minimum-hit-span INT  Require this span of exact anchors per target (0 disables)\n\
         --fallback-kmer INT  Query this shorter k only when earlier tiers have no hit; repeatable\n\
         --link-rad-arms      Route a hit to its sibling __R1/__R2 arm\n\
         --link-rad-max-fragments INT  Stop sibling propagation after this many direct fragments (default: 256)\n\
         --version     Print version"
    );
}

// 取紧跟在开关后头的值；少了就当场吱声，别往下带病跑。
fn take_value(argv: &[String], index: &mut usize, option: &str) -> AppResult<String> {
    *index += 1;
    argv.get(*index)
        .cloned()
        .ok_or_else(|| format!("option {option} requires an argument"))
}

// 数字参数统一在这儿验，免得各处转来转去整乱套。
fn parse_usize(value: String, option: &str) -> AppResult<usize> {
    value
        .parse::<usize>()
        .map_err(|_| format!("cannot convert {value:?} for {option} to an integer"))
}

fn parse_args(argv: Vec<String>) -> AppResult<Args> {
    // 先把参数归拢明白，后面干活儿就不抓瞎了。
    if argv.iter().any(|arg| arg == "-h" || arg == "--help") {
        print_help();
        process::exit(0);
    }
    if argv.iter().any(|arg| arg == "--version") {
        println!("MainFilterNew 0.4.0");
        process::exit(0);
    }

    let mut parsed = Args::default();
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-r" => parsed.reference = PathBuf::from(take_value(&argv, &mut i, "-r")?),
            "-o" => parsed.output = PathBuf::from(take_value(&argv, &mut i, "-o")?),
            "-kf" => parsed.kmer = parse_usize(take_value(&argv, &mut i, "-kf")?, "-kf")?,
            "-s" => parsed.step = parse_usize(take_value(&argv, &mut i, "-s")?, "-s")?,
            "-m_reads" => {
                parsed.max_read_blocks = take_value(&argv, &mut i, "-m_reads")?
                    .parse::<u64>()
                    .map_err(|_| "-m_reads requires a non-negative integer".to_string())?;
            }
            "-subdir" => parsed.out_subdir = take_value(&argv, &mut i, "-subdir")?,
            "-lkd" => parsed.dictionary = Some(PathBuf::from(take_value(&argv, &mut i, "-lkd")?)),
            "-m" => {
                parsed.mode = take_value(&argv, &mut i, "-m")?
                    .parse::<u8>()
                    .map_err(|_| "-m requires an integer from 0 to 5".to_string())?;
            }
            "--minimum-hit-span" => {
                parsed.minimum_hit_span = parse_usize(
                    take_value(&argv, &mut i, "--minimum-hit-span")?,
                    "--minimum-hit-span",
                )?
            }
            "--fallback-kmer" => parsed.fallback_kmers.push(parse_usize(
                take_value(&argv, &mut i, "--fallback-kmer")?,
                "--fallback-kmer",
            )?),
            "--link-rad-arms" => parsed.link_rad_arms = true,
            "--link-rad-max-fragments" => {
                parsed.link_rad_max_fragments = parse_usize(
                    take_value(&argv, &mut i, "--link-rad-max-fragments")?,
                    "--link-rad-max-fragments",
                )? as u64;
            }
            "-gr" => parsed.get_reverse = true,
            "-lb" => parsed.use_composition_pattern = true,
            "-q1" | "-q2" => {
                let is_q1 = argv[i] == "-q1";
                let mut paths = Vec::new();
                while i + 1 < argv.len() && !argv[i + 1].starts_with('-') {
                    i += 1;
                    paths.push(PathBuf::from(&argv[i]));
                }
                if paths.is_empty() {
                    return Err(format!(
                        "option {} requires at least one file",
                        if is_q1 { "-q1" } else { "-q2" }
                    ));
                }
                if is_q1 {
                    parsed.q1.extend(paths);
                } else {
                    parsed.q2.extend(paths);
                }
            }
            option if option.starts_with('-') => return Err(format!("unknown option {option}")),
            value => return Err(format!("unexpected argument {value}")),
        }
        i += 1;
    }

    if parsed.reference.as_os_str().is_empty() {
        return Err("reference sequences must be supplied".to_string());
    }
    if parsed.output.as_os_str().is_empty() {
        return Err("output directory must be supplied".to_string());
    }
    if parsed.kmer < 16 {
        return Err("k-mer size must be at least 16".to_string());
    }
    if parsed.step == 0 {
        return Err("read scanning step must be at least 1".to_string());
    }
    if parsed.mode > 5 {
        return Err("mode must be between 0 and 5".to_string());
    }
    let mut previous = parsed.kmer;
    for &fallback in &parsed.fallback_kmers {
        if fallback < 16 || fallback >= previous {
            return Err("fallback k-mers must be at least 16 and strictly decreasing".to_string());
        }
        previous = fallback;
    }
    if !parsed.q2.is_empty() && parsed.q1.len() != parsed.q2.len() {
        return Err("-q1 and -q2 must contain the same number of files".to_string());
    }
    if parsed.mode != 2 && parsed.q1.is_empty() {
        return Err("at least one sequencing file is required".to_string());
    }
    Ok(parsed)
}

// 系统 zlib 静态链进来,当兜底;zlib-ng 有没有得等运行时探测才知道。
#[link(name = "z")]
extern "C" {
    fn gzopen(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn gzread(file: *mut c_void, buffer: *mut c_void, length: u32) -> c_int;
    fn gzclose(file: *mut c_void) -> c_int;
    fn gzbuffer(file: *mut c_void, size: c_uint) -> c_int;
}

type GzOpenFn = unsafe extern "C" fn(*const c_char, *const c_char) -> *mut c_void;
type GzReadFn = unsafe extern "C" fn(*mut c_void, *mut c_void, u32) -> c_int;
type GzCloseFn = unsafe extern "C" fn(*mut c_void) -> c_int;
type GzBufferFn = unsafe extern "C" fn(*mut c_void, c_uint) -> c_int;

#[derive(Clone, Copy)]
struct ZlibBackend {
    open: GzOpenFn,
    read: GzReadFn,
    close: GzCloseFn,
    buffer: GzBufferFn,
    name: &'static str,
}

fn stock_zlib_backend() -> ZlibBackend {
    ZlibBackend {
        open: gzopen,
        read: gzread,
        close: gzclose,
        buffer: gzbuffer,
        name: "system zlib",
    }
}

// 构建时经 pkg-config 确认的原生 zlib-ng 直接链接，避免每次启动再走
// dlopen/dlsym。保留运行时探测作为可移植的后备路径。
#[cfg(system_zlib_ng)]
extern "C" {
    fn zng_gzopen(path: *const c_char, mode: *const c_char) -> *mut c_void;
    fn zng_gzread(file: *mut c_void, buffer: *mut c_void, length: u32) -> c_int;
    fn zng_gzclose(file: *mut c_void) -> c_int;
    fn zng_gzbuffer(file: *mut c_void, size: c_uint) -> c_int;
}

#[cfg(system_zlib_ng)]
fn build_detected_zlib_ng_backend() -> Option<ZlibBackend> {
    Some(ZlibBackend {
        open: zng_gzopen,
        read: zng_gzread,
        close: zng_gzclose,
        buffer: zng_gzbuffer,
        name: "zlib-ng (build detected)",
    })
}

#[cfg(not(system_zlib_ng))]
fn build_detected_zlib_ng_backend() -> Option<ZlibBackend> {
    None
}

// dlopen 一个符号,类型对不上就当没找到,不瞎猜。
#[cfg(unix)]
/// # Safety
///
/// `handle` 必须是仍然有效的 `dlopen` 句柄，且 `F` 必须与 `symbol` 的
/// 导出 ABI 和函数签名完全一致。调用方在使用返回的函数指针期间保持库已加载。
unsafe fn dlsym_typed<F: Copy>(handle: *mut c_void, symbol: &str) -> Option<F> {
    let symbol = CString::new(symbol).ok()?;
    let pointer = libc::dlsym(handle, symbol.as_ptr());
    if pointer.is_null() {
        None
    } else {
        Some(std::mem::transmute_copy(&pointer))
    }
}

// 运行时找 zlib-ng(原生 API,zng_ 前缀符号):装了就用它的 SIMD 加速解压,
// 没装、或者库里缺符号,就用上面静态链的系统 zlib,两者对上层完全透明。
#[cfg(unix)]
fn detect_zlib_ng() -> Option<ZlibBackend> {
    const CANDIDATE_LIBRARY_NAMES: &[&str] = &["libz-ng.so.2", "libz-ng.so.1", "libz-ng.so"];
    for library_name in CANDIDATE_LIBRARY_NAMES {
        let library_name = CString::new(*library_name).expect("static string contains no NUL");
        let handle =
            unsafe { libc::dlopen(library_name.as_ptr(), libc::RTLD_NOW | libc::RTLD_LOCAL) };
        if handle.is_null() {
            continue;
        }
        let resolved = unsafe {
            (
                dlsym_typed::<GzOpenFn>(handle, "zng_gzopen"),
                dlsym_typed::<GzReadFn>(handle, "zng_gzread"),
                dlsym_typed::<GzCloseFn>(handle, "zng_gzclose"),
                dlsym_typed::<GzBufferFn>(handle, "zng_gzbuffer"),
            )
        };
        if let (Some(open), Some(read), Some(close), Some(buffer)) = resolved {
            // 故意不 dlclose:句柄要活到进程退出,函数指针才一直有效。
            return Some(ZlibBackend {
                open,
                read,
                close,
                buffer,
                name: "zlib-ng",
            });
        }
        unsafe { libc::dlclose(handle) };
    }
    None
}

#[cfg(not(unix))]
fn detect_zlib_ng() -> Option<ZlibBackend> {
    None
}

static ZLIB_BACKEND: std::sync::OnceLock<ZlibBackend> = std::sync::OnceLock::new();

fn zlib_backend() -> ZlibBackend {
    *ZLIB_BACKEND.get_or_init(|| {
        build_detected_zlib_ng_backend()
            .or_else(detect_zlib_ng)
            .unwrap_or_else(stock_zlib_backend)
    })
}

struct GzipReader {
    handle: *mut c_void,
    backend: ZlibBackend,
}

impl GzipReader {
    fn open(path: &Path) -> io::Result<Self> {
        let backend = zlib_backend();
        let path = CString::new(path.to_string_lossy().as_bytes())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains a NUL byte"))?;
        let mode = CString::new("rb").expect("static string contains no NUL");
        let handle = unsafe { (backend.open)(path.as_ptr(), mode.as_ptr()) };
        if handle.is_null() {
            Err(io::Error::new(
                io::ErrorKind::NotFound,
                "cannot open gzip file",
            ))
        } else {
            // 把压缩端缓冲从默认 8KiB 拉大,少喊几次底层 read(2)。
            // 必须在第一次 gzread 之前调用,失败也无妨,大不了退回默认缓冲区大小。
            unsafe { (backend.buffer)(handle, READ_BUFFER_SIZE as c_uint) };
            Ok(Self { handle, backend })
        }
    }
}

impl Read for GzipReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        let size = buffer.len().min(c_int::MAX as usize) as u32;
        let result = unsafe { (self.backend.read)(self.handle, buffer.as_mut_ptr().cast(), size) };
        if result < 0 {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "gzip decompression failed",
            ))
        } else {
            Ok(result as usize)
        }
    }
}

impl Drop for GzipReader {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            unsafe { (self.backend.close)(self.handle) };
        }
    }
}

// 看文件头认 gzip，比光瞅后缀靠谱点儿。
fn is_gzip(path: &Path) -> bool {
    path.extension()
        .and_then(|value| value.to_str())
        .is_some_and(|value| value.eq_ignore_ascii_case("gz"))
}

fn sequence_extension(path: &Path) -> Option<String> {
    let base = if is_gzip(path) {
        path.file_stem().map(PathBuf::from)?
    } else {
        path.to_path_buf()
    };
    base.extension()
        .and_then(|value| value.to_str())
        .map(|value| value.to_ascii_lowercase())
}

// 先认清是 FASTA 还是 FASTQ，后面读法才不能串台。
fn detect_kind(path: &Path) -> AppResult<FileKind> {
    match sequence_extension(path).as_deref() {
        Some("fa" | "fas" | "fasta") => Ok(FileKind::Fasta),
        Some("fq" | "fastq") => Ok(FileKind::Fastq),
        _ if is_gzip(path) => Ok(FileKind::Fastq),
        _ => Err(format!(
            "unrecognized sequence file type: {}",
            path.display()
        )),
    }
}

// gzip 和普通文件都包成同一种读口，后头不用分两套流程。
fn open_input(path: &Path) -> io::Result<BufReader<Box<dyn Read>>> {
    let input: Box<dyn Read> = if is_gzip(path) {
        Box::new(GzipReader::open(path)?)
    } else {
        Box::new(File::open(path)?)
    };
    Ok(BufReader::with_capacity(READ_BUFFER_SIZE, input))
}

#[derive(Clone, Debug, Default)]
struct Record {
    header: Vec<u8>,
    sequence: Vec<u8>,
    plus: Vec<u8>,
    quality: Vec<u8>,
    has_quality: bool,
}

impl Record {
    fn reset(&mut self) {
        self.header.clear();
        self.sequence.clear();
        self.plus.clear();
        self.quality.clear();
        self.has_quality = false;
    }
}

struct SequenceReader {
    input: BufReader<Box<dyn Read>>,
    kind: FileKind,
    pending_header: Option<Vec<u8>>,
    line_scratch: Vec<u8>,
    header_scratch: Vec<u8>,
    plus_scratch: Vec<u8>,
    finished: bool,
    keep_text_lines: bool,
}

impl SequenceReader {
    fn open(path: &Path, kind: FileKind, keep_text_lines: bool) -> AppResult<Self> {
        Ok(Self {
            input: open_input(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?,
            kind,
            pending_header: None,
            line_scratch: Vec::with_capacity(512),
            header_scratch: Vec::with_capacity(256),
            plus_scratch: Vec::with_capacity(16),
            finished: false,
            keep_text_lines,
        })
    }

    fn read_line_into(&mut self, line: &mut Vec<u8>) -> AppResult<bool> {
        line.clear();
        if self
            .input
            .read_until(b'\n', line)
            .map_err(|e| e.to_string())?
            == 0
        {
            return Ok(false);
        }
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        Ok(true)
    }

    fn next_record_into(&mut self, record: &mut Record) -> AppResult<bool> {
        record.reset();
        match self.kind {
            FileKind::Fasta => self.next_fasta_into(record),
            FileKind::Fastq => self.next_fastq_into(record),
        }
    }

    fn read_scratch(&mut self) -> AppResult<bool> {
        let mut line = std::mem::take(&mut self.line_scratch);
        let present = self.read_line_into(&mut line)?;
        self.line_scratch = line;
        Ok(present)
    }

    fn next_fasta_into(&mut self, record: &mut Record) -> AppResult<bool> {
        if self.finished {
            return Ok(false);
        }
        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                if !self.read_scratch()? {
                    self.finished = true;
                    return Ok(false);
                }
                if self.line_scratch.first() == Some(&b'>') {
                    break std::mem::take(&mut self.line_scratch);
                }
                if !self.line_scratch.is_empty() {
                    return Err("FASTA sequence encountered before a header".to_string());
                }
            }
        };
        loop {
            if !self.read_scratch()? {
                self.finished = true;
                break;
            }
            if self.line_scratch.first() == Some(&b'>') {
                self.pending_header = Some(std::mem::take(&mut self.line_scratch));
                break;
            }
            record.sequence.extend_from_slice(&self.line_scratch);
        }
        if self.keep_text_lines {
            record.header = header;
        }
        Ok(true)
    }

    fn next_fastq_into(&mut self, record: &mut Record) -> AppResult<bool> {
        loop {
            let mut header = std::mem::take(&mut self.header_scratch);
            let present = self.read_line_into(&mut header)?;
            self.header_scratch = header;
            if !present {
                return Ok(false);
            }
            if !self.header_scratch.is_empty() {
                break;
            }
        }
        if self.header_scratch.first() != Some(&b'@') {
            return Err("malformed FASTQ record".to_string());
        }
        if self.keep_text_lines {
            record.header.extend_from_slice(&self.header_scratch);
        }
        if !self.read_line_into(&mut record.sequence)? {
            return Err("truncated FASTQ sequence".to_string());
        }
        let mut plus = std::mem::take(&mut self.plus_scratch);
        if !self.read_line_into(&mut plus)? {
            return Err("truncated FASTQ plus line".to_string());
        }
        self.plus_scratch = plus;
        if self.plus_scratch.first() != Some(&b'+') {
            return Err("malformed FASTQ plus line".to_string());
        }
        if self.keep_text_lines {
            record.plus.extend_from_slice(&self.plus_scratch);
        }
        if !self.read_line_into(&mut record.quality)? {
            return Err("truncated FASTQ quality".to_string());
        }
        if record.quality.len() != record.sequence.len() {
            return Err("FASTQ sequence and quality lengths differ".to_string());
        }
        record.has_quality = true;
        Ok(true)
    }
}

const INVALID_BASE_CODE: u8 = u8::MAX;

const fn build_base_code_table() -> [u8; 256] {
    let mut table = [INVALID_BASE_CODE; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
}

const BASE_CODE_TABLE: [u8; 256] = build_base_code_table();

#[inline(always)]
fn base_code(base: u8) -> Option<u8> {
    match BASE_CODE_TABLE[base as usize] {
        INVALID_BASE_CODE => None,
        code => Some(code),
    }
}

#[derive(Clone, Copy, Debug)]
enum ReferenceHits {
    One(u32),
    // Used only while building: each k-mer accumulates its locus IDs once.
    Pending(u32),
    // Final/read-only representation: a slice in packed_hits.
    Packed { offset: u32, len: u32 },
}

#[derive(Clone, Debug)]
enum KmerStore {
    Short(AHashMap<u64, ReferenceHits>),
    Medium(AHashMap<u128, ReferenceHits>),
    Long(AHashMap<Vec<u8>, ReferenceHits>),
}

impl KmerStore {
    fn new(k: usize) -> Self {
        match k {
            0..=32 => Self::Short(AHashMap::new()),
            33..=64 => Self::Medium(AHashMap::new()),
            _ => Self::Long(AHashMap::new()),
        }
    }
    fn len(&self) -> usize {
        match self {
            Self::Short(map) => map.len(),
            Self::Medium(map) => map.len(),
            Self::Long(map) => map.len(),
        }
    }
    fn reserve(&mut self, additional: usize) {
        match self {
            Self::Short(map) => map.reserve(additional),
            Self::Medium(map) => map.reserve(additional),
            Self::Long(map) => map.reserve(additional),
        }
    }
}

#[derive(Clone, Debug)]
struct KmerIndex {
    k: usize,
    reference_names: Vec<String>,
    store: KmerStore,
    pending_hits: Vec<Vec<u32>>,
    packed_hits: Vec<u32>,
}

impl KmerIndex {
    fn new(k: usize, reference_names: Vec<String>) -> Self {
        Self {
            k,
            reference_names,
            store: KmerStore::new(k),
            pending_hits: Vec::new(),
            packed_hits: Vec::new(),
        }
    }
    fn len(&self) -> usize {
        self.store.len()
    }

    fn add_reference_hit(
        hit: &mut ReferenceHits,
        reference: u32,
        pending_hits: &mut Vec<Vec<u32>>,
    ) {
        match *hit {
            ReferenceHits::One(previous) if previous != reference => {
                let index = pending_hits.len() as u32;
                pending_hits.push(vec![previous, reference]);
                *hit = ReferenceHits::Pending(index);
            }
            ReferenceHits::Pending(index) => {
                let hits = &mut pending_hits[index as usize];
                if hits.last().copied() != Some(reference) {
                    hits.push(reference);
                }
            }
            ReferenceHits::One(_) => {}
            ReferenceHits::Packed { .. } => unreachable!("cannot extend finalized reference hits"),
        }
    }

    fn insert_short(&mut self, key: u64, reference: u32) {
        let KmerStore::Short(map) = &mut self.store else {
            unreachable!("short k-mer inserted into long-kmer index");
        };
        match map.entry(key) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(ReferenceHits::One(reference));
            }
            std::collections::hash_map::Entry::Occupied(mut entry) => {
                Self::add_reference_hit(entry.get_mut(), reference, &mut self.pending_hits)
            }
        }
    }

    fn insert_medium(&mut self, key: u128, reference: u32) {
        let KmerStore::Medium(map) = &mut self.store else {
            unreachable!("medium k-mer inserted into a different index");
        };
        match map.entry(key) {
            std::collections::hash_map::Entry::Vacant(entry) => {
                entry.insert(ReferenceHits::One(reference));
            }
            std::collections::hash_map::Entry::Occupied(mut entry) => {
                Self::add_reference_hit(entry.get_mut(), reference, &mut self.pending_hits)
            }
        }
    }

    fn insert_long(&mut self, key: &[u8], reference: u32) {
        let KmerStore::Long(map) = &mut self.store else {
            unreachable!("long k-mer inserted into short-kmer index");
        };
        if let Some(hit) = map.get_mut(key) {
            Self::add_reference_hit(hit, reference, &mut self.pending_hits);
        } else {
            map.insert(key.to_vec(), ReferenceHits::One(reference));
        }
    }

    fn pack_pending_hit(
        hit: &mut ReferenceHits,
        pending_hits: &[Vec<u32>],
        packed_hits: &mut Vec<u32>,
    ) -> AppResult<()> {
        if let ReferenceHits::Pending(index) = *hit {
            let values = &pending_hits[index as usize];
            let offset = u32::try_from(packed_hits.len())
                .map_err(|_| "packed reference hits exceed u32 range".to_string())?;
            let len = u32::try_from(values.len())
                .map_err(|_| "reference hit list exceeds u32 range".to_string())?;
            offset
                .checked_add(len)
                .ok_or_else(|| "packed reference hit range exceeds u32 range".to_string())?;
            packed_hits.extend_from_slice(values);
            *hit = ReferenceHits::Packed { offset, len };
        }
        Ok(())
    }

    fn finalize_hits(&mut self) -> AppResult<()> {
        match &mut self.store {
            KmerStore::Short(map) => {
                for hit in map.values_mut() {
                    Self::pack_pending_hit(hit, &self.pending_hits, &mut self.packed_hits)?;
                }
            }
            KmerStore::Medium(map) => {
                for hit in map.values_mut() {
                    Self::pack_pending_hit(hit, &self.pending_hits, &mut self.packed_hits)?;
                }
            }
            KmerStore::Long(map) => {
                for hit in map.values_mut() {
                    Self::pack_pending_hit(hit, &self.pending_hits, &mut self.packed_hits)?;
                }
            }
        }
        self.pending_hits.clear();
        self.pending_hits.shrink_to_fit();
        Ok(())
    }

    fn hit_slice<'a>(&'a self, hits: &'a ReferenceHits) -> &'a [u32] {
        match hits {
            ReferenceHits::One(reference) => std::slice::from_ref(reference),
            ReferenceHits::Packed { offset, len } => {
                let end = offset
                    .checked_add(*len)
                    .expect("validated packed reference-hit range");
                &self.packed_hits[*offset as usize..end as usize]
            }
            ReferenceHits::Pending(_) => {
                unreachable!("reference hits must be finalized before querying")
            }
        }
    }

    fn short_mask(&self) -> u64 {
        if self.k == 32 {
            u64::MAX
        } else {
            (1_u64 << (2 * self.k)) - 1
        }
    }

    fn medium_mask(&self) -> u128 {
        if self.k == 64 {
            u128::MAX
        } else {
            (1_u128 << (2 * self.k)) - 1
        }
    }

    fn add_reference_sequence(&mut self, sequence: &[u8], reference: u32) {
        if sequence.len() < self.k {
            return;
        }
        match self.k {
            0..=32 => {
                let mask = self.short_mask();
                let reverse_shift = 2 * (self.k - 1);
                let (mut forward, mut reverse, mut valid) = (0_u64, 0_u64, 0_usize);
                for &base in sequence {
                    if let Some(code) = base_code(base) {
                        forward = ((forward << 2) | code as u64) & mask;
                        reverse = (reverse >> 2) | (((3 - code) as u64) << reverse_shift);
                        valid += 1;
                        if valid >= self.k {
                            self.insert_short(forward.min(reverse), reference);
                        }
                    } else {
                        forward = 0;
                        reverse = 0;
                        valid = 0;
                    }
                }
            }
            33..=64 => {
                let mask = self.medium_mask();
                let reverse_shift = 2 * (self.k - 1);
                let (mut forward, mut reverse, mut valid) = (0_u128, 0_u128, 0_usize);
                for &base in sequence {
                    if let Some(code) = base_code(base) {
                        forward = ((forward << 2) | code as u128) & mask;
                        reverse = (reverse >> 2) | (((3 - code) as u128) << reverse_shift);
                        valid += 1;
                        if valid >= self.k {
                            self.insert_medium(forward.min(reverse), reference);
                        }
                    } else {
                        forward = 0;
                        reverse = 0;
                        valid = 0;
                    }
                }
            }
            _ => {
                let mut forward = Vec::with_capacity(self.k);
                let mut reverse = Vec::with_capacity(self.k);
                for start in 0..=sequence.len() - self.k {
                    if long_kmer_into(&sequence[start..start + self.k], &mut forward, &mut reverse)
                    {
                        self.insert_long(
                            if forward <= reverse {
                                &forward
                            } else {
                                &reverse
                            },
                            reference,
                        );
                    }
                }
            }
        }
    }

    fn collect_hits(&self, sequence: &[u8], step: usize, collector: &mut HitCollector) {
        if sequence.len() < self.k {
            return;
        }
        let tail = sequence.len() - self.k;
        match self.k {
            0..=32 => {
                let mask = self.short_mask();
                let reverse_shift = 2 * (self.k - 1);
                let (mut forward, mut reverse, mut valid, mut next_probe) =
                    (0_u64, 0_u64, 0_usize, 0_usize);
                for (end, &base) in sequence.iter().enumerate() {
                    if let Some(code) = base_code(base) {
                        forward = ((forward << 2) | code as u64) & mask;
                        reverse = (reverse >> 2) | (((3 - code) as u64) << reverse_shift);
                        valid += 1;
                    } else {
                        forward = 0;
                        reverse = 0;
                        valid = 0;
                    }
                    if end + 1 < self.k {
                        continue;
                    }
                    let start = end + 1 - self.k;
                    let sampled = start == next_probe;
                    if sampled {
                        next_probe = next_probe.saturating_add(step);
                    }
                    if valid >= self.k && (sampled || start == tail) {
                        collector.position = collector.origin + start;
                        self.collect_short(forward.min(reverse), collector);
                    }
                }
            }
            33..=64 => {
                let mask = self.medium_mask();
                let reverse_shift = 2 * (self.k - 1);
                let (mut forward, mut reverse, mut valid, mut next_probe) =
                    (0_u128, 0_u128, 0_usize, 0_usize);
                for (end, &base) in sequence.iter().enumerate() {
                    if let Some(code) = base_code(base) {
                        forward = ((forward << 2) | code as u128) & mask;
                        reverse = (reverse >> 2) | (((3 - code) as u128) << reverse_shift);
                        valid += 1;
                    } else {
                        forward = 0;
                        reverse = 0;
                        valid = 0;
                    }
                    if end + 1 < self.k {
                        continue;
                    }
                    let start = end + 1 - self.k;
                    let sampled = start == next_probe;
                    if sampled {
                        next_probe = next_probe.saturating_add(step);
                    }
                    if valid >= self.k && (sampled || start == tail) {
                        collector.position = collector.origin + start;
                        self.collect_medium(forward.min(reverse), collector);
                    }
                }
            }
            _ => {
                let mut forward = Vec::with_capacity(self.k);
                let mut reverse = Vec::with_capacity(self.k);
                for start in (0..=tail)
                    .step_by(step)
                    .chain((!tail.is_multiple_of(step)).then_some(tail))
                {
                    if long_kmer_into(&sequence[start..start + self.k], &mut forward, &mut reverse)
                    {
                        collector.position = collector.origin + start;
                        self.collect_long(
                            if forward <= reverse {
                                &forward
                            } else {
                                &reverse
                            },
                            collector,
                        );
                    }
                }
            }
        }
    }

    fn collect_short(&self, key: u64, collector: &mut HitCollector) {
        let KmerStore::Short(map) = &self.store else {
            unreachable!("short lookup in long-kmer index");
        };
        if let Some(hits) = map.get(&key) {
            for &reference in self.hit_slice(hits) {
                collector.mark(reference as usize);
            }
        }
    }

    fn collect_medium(&self, key: u128, collector: &mut HitCollector) {
        let KmerStore::Medium(map) = &self.store else {
            unreachable!("medium lookup in a different index");
        };
        if let Some(hits) = map.get(&key) {
            for &reference in self.hit_slice(hits) {
                collector.mark(reference as usize);
            }
        }
    }

    fn collect_long(&self, key: &[u8], collector: &mut HitCollector) {
        let KmerStore::Long(map) = &self.store else {
            unreachable!("long lookup in short-kmer index");
        };
        if let Some(hits) = map.get(key) {
            for &reference in self.hit_slice(hits) {
                collector.mark(reference as usize);
            }
        }
    }
}

// k 太长塞不进 u64 时，老老实实造正反链字符串键。
fn long_kmer_into(sequence: &[u8], forward: &mut Vec<u8>, reverse: &mut Vec<u8>) -> bool {
    forward.clear();
    reverse.clear();
    if forward.capacity() < sequence.len() {
        forward.reserve(sequence.len() - forward.capacity());
    }
    if reverse.capacity() < sequence.len() {
        reverse.reserve(sequence.len() - reverse.capacity());
    }
    for &base in sequence {
        let Some(code) = base_code(base) else {
            forward.clear();
            reverse.clear();
            return false;
        };
        forward.push(b"ACGT"[code as usize]);
    }
    for &base in sequence.iter().rev() {
        let Some(code) = base_code(base) else {
            unreachable!("validated above");
        };
        reverse.push(b"ACGT"[(3 - code) as usize]);
    }
    true
}

struct HitCollector {
    generation: u32,
    seen: Vec<u32>,
    hits: Vec<u32>,
    first: Vec<usize>,
    last: Vec<usize>,
    position: usize,
    origin: usize,
    minimum_span: usize,
    required_span: usize,
}

impl HitCollector {
    fn new(reference_count: usize) -> Self {
        Self {
            generation: 0,
            seen: vec![0; reference_count],
            hits: Vec::new(),
            first: vec![0; reference_count],
            last: vec![0; reference_count],
            position: 0,
            origin: 0,
            minimum_span: 0,
            required_span: 0,
        }
    }

    fn begin(&mut self) {
        // 换个批次号就当清空了，省得每条 read 都把整张表擦一遍。
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.seen.fill(0);
            self.generation = 1;
        }
        self.hits.clear();
    }

    fn mark(&mut self, reference: usize) {
        if self.seen[reference] != self.generation {
            self.seen[reference] = self.generation;
            self.hits.push(reference as u32);
            if self.required_span > 0 {
                self.first[reference] = self.position;
                self.last[reference] = self.position;
            }
        } else if self.required_span > 0
            && self.last[reference] - self.first[reference] < self.required_span
        {
            self.first[reference] = self.first[reference].min(self.position);
            self.last[reference] = self.last[reference].max(self.position);
        }
    }

    fn link_siblings(&mut self, siblings: &[Option<u32>], direct_counts: &mut [u64], maximum: u64) {
        let direct_hits = self.hits.len();
        for index in 0..direct_hits {
            let reference = self.hits[index] as usize;
            direct_counts[reference] = direct_counts[reference].saturating_add(1);
            if direct_counts[reference] <= maximum {
                if let Some(sibling) = siblings[reference] {
                    self.mark(sibling as usize);
                }
            }
        }
    }
}

fn rad_arm_siblings(names: &[String]) -> AppResult<Vec<Option<u32>>> {
    let mut arms = AHashMap::<&str, [Option<usize>; 2]>::new();
    for (index, name) in names.iter().enumerate() {
        let (locus, arm) = if let Some(locus) = name.strip_suffix("__R1") {
            (locus, 0)
        } else if let Some(locus) = name.strip_suffix("__R2") {
            (locus, 1)
        } else {
            continue;
        };
        let slot = &mut arms.entry(locus).or_default()[arm];
        if slot.replace(index).is_some() {
            return Err(format!("duplicate RAD arm reference: {name}"));
        }
    }
    let mut siblings = vec![None; names.len()];
    for pair in arms.values() {
        if let [Some(left), Some(right)] = *pair {
            siblings[left] = Some(right as u32);
            siblings[right] = Some(left as u32);
        }
    }
    Ok(siblings)
}

fn collect_fragment_hits(
    primary: &KmerIndex,
    fallbacks: &[KmerIndex],
    first: &[u8],
    second: Option<&[u8]>,
    step: usize,
    collector: &mut HitCollector,
) {
    for index in std::iter::once(primary).chain(fallbacks.iter()) {
        collector.required_span = collector.minimum_span.saturating_sub(index.k);
        collector.origin = 0;
        index.collect_hits(first, step, collector);
        if let Some(second) = second {
            // Separate mates: one anchor on each is independent evidence.
            collector.origin = first.len() + index.k;
            index.collect_hits(second, step, collector);
        }
        if collector.minimum_span > index.k {
            collector.hits.retain(|&id| {
                collector.last[id as usize] - collector.first[id as usize] + index.k
                    >= collector.minimum_span
            });
        }
        if !collector.hits.is_empty() {
            break;
        }
        // A rejected short anchor cannot leak support into another tier.
        collector.begin();
    }
}

// 文件名就是 locus 名，去掉 .gz 和序列后缀，别把扩展名带进结果。
fn reference_basename(path: &Path) -> AppResult<String> {
    let without_gzip = if is_gzip(path) {
        path.file_stem()
            .map(PathBuf::from)
            .ok_or_else(|| "invalid reference name".to_string())?
    } else {
        path.to_path_buf()
    };
    without_gzip
        .file_stem()
        .and_then(|value| value.to_str())
        .map(str::to_string)
        .ok_or_else(|| format!("invalid reference name: {}", path.display()))
}

// 参考既能给一个文件，也能给一摞文件；这儿统一排好队。
fn reference_paths(reference: &Path) -> AppResult<Vec<PathBuf>> {
    let mut paths = Vec::new();
    if reference.is_dir() {
        for entry in fs::read_dir(reference).map_err(|e| e.to_string())? {
            let path = entry.map_err(|e| e.to_string())?.path();
            if path.is_file() && detect_kind(&path).ok() == Some(FileKind::Fasta) {
                paths.push(path);
            }
        }
    } else if reference.is_file() {
        paths.push(reference.to_path_buf());
    }
    paths.sort();
    if paths.is_empty() {
        return Err("no reference FASTA file found".to_string());
    }
    Ok(paths)
}

fn reference_content_hash(reference: &Path) -> AppResult<[u8; 32]> {
    let paths = reference_paths(reference)?;
    let mut hasher = Sha256::new();
    hasher.update(b"GM2-MainFilter-reference-v1\0");
    for path in paths {
        let name = reference_basename(&path)?;
        hasher.update((name.len() as u64).to_le_bytes());
        hasher.update(name.as_bytes());
        let mut reader = SequenceReader::open(&path, FileKind::Fasta, false)?;
        let mut record = Record::default();
        while reader.next_record_into(&mut record)? {
            hasher.update((record.sequence.len() as u64).to_le_bytes());
            hasher.update(&record.sequence);
        }
    }
    Ok(hasher.finalize().into())
}

fn build_index(
    reference: &Path,
    k: usize,
    compute_hash: bool,
) -> AppResult<(KmerIndex, Option<[u8; 32]>)> {
    // 首次构建时把内容哈希和索引合在一遍参考读取里，避免额外 I/O。
    let paths = reference_paths(reference)?;
    let names: Vec<String> = paths
        .iter()
        .map(|path| reference_basename(path))
        .collect::<AppResult<_>>()?;
    let mut unique = HashSet::new();
    for name in &names {
        if !unique.insert(name.clone()) {
            return Err(format!("duplicate reference locus name: {name}"));
        }
    }
    let mut hasher = compute_hash.then(Sha256::new);
    if let Some(hasher) = &mut hasher {
        hasher.update(b"GM2-MainFilter-reference-v1\0");
    }
    let mut index = KmerIndex::new(k, names);
    for (reference_id, path) in paths.iter().enumerate() {
        let name = &index.reference_names[reference_id];
        if let Some(hasher) = &mut hasher {
            hasher.update((name.len() as u64).to_le_bytes());
            hasher.update(name.as_bytes());
        }
        let mut reader = SequenceReader::open(path, FileKind::Fasta, false)?;
        let mut record = Record::default();
        while reader.next_record_into(&mut record)? {
            if let Some(hasher) = &mut hasher {
                hasher.update((record.sequence.len() as u64).to_le_bytes());
                hasher.update(&record.sequence);
            }
            index.add_reference_sequence(&record.sequence, reference_id as u32);
        }
    }
    index.finalize_hits()?;
    Ok((index, hasher.map(|hasher| hasher.finalize().into())))
}

fn write_u16(out: &mut impl Write, value: u16) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}
fn write_u32(out: &mut impl Write, value: u32) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}
fn write_u64(out: &mut impl Write, value: u64) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}
fn write_u128(out: &mut impl Write, value: u128) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}

fn write_dictionary(index: &KmerIndex, reference_hash: &[u8; 32], path: &Path) -> AppResult<()> {
    // 先写临时文件再改名，半道断了也不至于把老缓存整坏喽。
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let temporary = path.with_extension(format!("tmp.{}", process::id()));
    let file = File::create(&temporary).map_err(|e| e.to_string())?;
    let mut out = BufWriter::with_capacity(4 * 1024 * 1024, file);
    out.write_all(CACHE_MAGIC).map_err(|e| e.to_string())?;
    write_u16(&mut out, CACHE_VERSION).map_err(|e| e.to_string())?;
    write_u16(&mut out, FLAG_CANONICAL_KMERS).map_err(|e| e.to_string())?;
    write_u32(&mut out, index.k as u32).map_err(|e| e.to_string())?;
    write_u32(&mut out, index.reference_names.len() as u32).map_err(|e| e.to_string())?;
    write_u64(&mut out, index.len() as u64).map_err(|e| e.to_string())?;
    out.write_all(reference_hash).map_err(|e| e.to_string())?;
    for name in &index.reference_names {
        write_u32(&mut out, name.len() as u32).map_err(|e| e.to_string())?;
        out.write_all(name.as_bytes()).map_err(|e| e.to_string())?;
    }
    match &index.store {
        KmerStore::Short(map) => {
            for (value, hits) in map {
                out.write_all(&[0]).map_err(|e| e.to_string())?;
                write_u64(&mut out, *value).map_err(|e| e.to_string())?;
                let hits = index.hit_slice(hits);
                write_u32(&mut out, hits.len() as u32).map_err(|e| e.to_string())?;
                for &hit in hits {
                    write_u32(&mut out, hit).map_err(|e| e.to_string())?;
                }
            }
        }
        KmerStore::Medium(map) => {
            for (value, hits) in map {
                out.write_all(&[1]).map_err(|e| e.to_string())?;
                write_u128(&mut out, *value).map_err(|e| e.to_string())?;
                let hits = index.hit_slice(hits);
                write_u32(&mut out, hits.len() as u32).map_err(|e| e.to_string())?;
                for &hit in hits {
                    write_u32(&mut out, hit).map_err(|e| e.to_string())?;
                }
            }
        }
        KmerStore::Long(map) => {
            for (value, hits) in map {
                out.write_all(&[2]).map_err(|e| e.to_string())?;
                write_u32(&mut out, value.len() as u32).map_err(|e| e.to_string())?;
                out.write_all(value).map_err(|e| e.to_string())?;
                let hits = index.hit_slice(hits);
                write_u32(&mut out, hits.len() as u32).map_err(|e| e.to_string())?;
                for &hit in hits {
                    write_u32(&mut out, hit).map_err(|e| e.to_string())?;
                }
            }
        }
    }
    out.flush().map_err(|e| e.to_string())?;
    fs::rename(&temporary, path).map_err(|e| e.to_string())?;
    Ok(())
}

fn read_array<const N: usize>(input: &mut impl Read) -> AppResult<[u8; N]> {
    let mut bytes = [0_u8; N];
    input
        .read_exact(&mut bytes)
        .map_err(|_| "truncated k-mer dictionary".to_string())?;
    Ok(bytes)
}

fn read_u8(input: &mut impl Read) -> AppResult<u8> {
    Ok(read_array::<1>(input)?[0])
}

fn read_u16(input: &mut impl Read) -> AppResult<u16> {
    Ok(u16::from_le_bytes(read_array(input)?))
}

fn read_u32(input: &mut impl Read) -> AppResult<u32> {
    Ok(u32::from_le_bytes(read_array(input)?))
}

fn read_u64(input: &mut impl Read) -> AppResult<u64> {
    Ok(u64::from_le_bytes(read_array(input)?))
}
fn read_u128(input: &mut impl Read) -> AppResult<u128> {
    Ok(u128::from_le_bytes(read_array(input)?))
}

fn read_vec(input: &mut impl Read, length: usize) -> AppResult<Vec<u8>> {
    let mut bytes = vec![0_u8; length];
    input
        .read_exact(&mut bytes)
        .map_err(|_| "truncated k-mer dictionary".to_string())?;
    Ok(bytes)
}

// 先看固定大小的头；版本或 k 不对时无需扫描整个参考计算内容哈希。
fn probe_dictionary(path: &Path, requested_k: usize) -> AppResult<()> {
    let file = File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut input = BufReader::with_capacity(64 * 1024, file);
    if &read_array::<4>(&mut input)? != CACHE_MAGIC {
        return Err("legacy or invalid k-mer dictionary; rebuilding".to_string());
    }
    if read_u16(&mut input)? != CACHE_VERSION {
        return Err("unsupported k-mer dictionary version".to_string());
    }
    let flags = read_u16(&mut input)?;
    if flags & FLAG_CANONICAL_KMERS == 0 {
        return Err("dictionary does not use canonical k-mers".to_string());
    }
    let stored_k = read_u32(&mut input)? as usize;
    if stored_k != requested_k {
        return Err(format!(
            "requested {requested_k}-mer but dictionary contains {stored_k}-mer"
        ));
    }
    Ok(())
}

// 优先尝试捡现成字典；格式和参数对不上就拒绝，不能瞎凑合。
fn load_dictionary(
    path: &Path,
    requested_k: usize,
    expected_reference_hash: &[u8; 32],
) -> AppResult<KmerIndex> {
    let file = File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut input = BufReader::new(file);
    let prefix = read_array::<4>(&mut input)?;
    if &prefix != CACHE_MAGIC {
        return Err(
            "legacy k-mer dictionary lacks a reference-content hash; rebuilding".to_string(),
        );
    }
    load_v3_dictionary(&mut input, requested_k, expected_reference_hash)
}

// 读取 Rust 缓存，核对 reference 内容、k 与 canonical 链策略，避免静默拿错字典。
fn load_v3_dictionary(
    input: &mut impl BufRead,
    requested_k: usize,
    expected_reference_hash: &[u8; 32],
) -> AppResult<KmerIndex> {
    if read_u16(input)? != CACHE_VERSION {
        return Err("unsupported k-mer dictionary version".to_string());
    }
    let flags = read_u16(input)?;
    if flags & FLAG_CANONICAL_KMERS == 0 {
        return Err("dictionary does not use canonical k-mers".to_string());
    }
    let k = read_u32(input)? as usize;
    if k != requested_k {
        return Err(format!(
            "requested {requested_k}-mer but dictionary contains {k}-mer"
        ));
    }
    let reference_count = read_u32(input)? as usize;
    let entry_count = usize::try_from(read_u64(input)?)
        .map_err(|_| "dictionary entry count exceeds this platform".to_string())?;
    let found_reference_hash = read_array::<32>(input)?;
    if &found_reference_hash != expected_reference_hash {
        return Err("dictionary reference-content hash does not match; rebuilding".to_string());
    }
    let mut names = Vec::with_capacity(reference_count);
    for _ in 0..reference_count {
        let length = read_u32(input)? as usize;
        names.push(
            String::from_utf8(read_vec(input, length)?)
                .map_err(|_| "reference name is not UTF-8".to_string())?,
        );
    }
    let mut index = KmerIndex::new(k, names);
    index.store.reserve(entry_count);
    for _ in 0..entry_count {
        let key_type = read_u8(input)?;
        match key_type {
            0 if k <= 32 => {
                let key = read_u64(input)?;
                let value =
                    read_dictionary_hits_into(input, reference_count, &mut index.packed_hits)?;
                let KmerStore::Short(map) = &mut index.store else {
                    unreachable!()
                };
                if map.insert(key, value).is_some() {
                    return Err("duplicate k-mer in dictionary".to_string());
                }
            }
            1 if (33..=64).contains(&k) => {
                let key = read_u128(input)?;
                let value =
                    read_dictionary_hits_into(input, reference_count, &mut index.packed_hits)?;
                let KmerStore::Medium(map) = &mut index.store else {
                    unreachable!()
                };
                if map.insert(key, value).is_some() {
                    return Err("duplicate k-mer in dictionary".to_string());
                }
            }
            2 if k > 64 => {
                let length = read_u32(input)? as usize;
                if length != k {
                    return Err("dictionary k-mer type or length does not match k".to_string());
                }
                let key = read_vec(input, length)?;
                let value =
                    read_dictionary_hits_into(input, reference_count, &mut index.packed_hits)?;
                let KmerStore::Long(map) = &mut index.store else {
                    unreachable!()
                };
                if map.insert(key, value).is_some() {
                    return Err("duplicate k-mer in dictionary".to_string());
                }
            }
            _ => return Err("dictionary k-mer type does not match k".to_string()),
        }
    }
    let mut trailing = [0_u8; 1];
    if input.read(&mut trailing).map_err(|e| e.to_string())? != 0 {
        return Err("unexpected trailing dictionary data".to_string());
    }
    Ok(index)
}

// 直接将多命中写入最终 packed_hits，避免每个 k-mer 分配临时 Vec。
fn read_dictionary_hits_into(
    input: &mut impl Read,
    reference_count: usize,
    packed_hits: &mut Vec<u32>,
) -> AppResult<ReferenceHits> {
    let hit_count = read_u32(input)? as usize;
    if hit_count == 0 {
        return Err("dictionary k-mer has no reference hits".to_string());
    }
    let offset = u32::try_from(packed_hits.len())
        .map_err(|_| "packed reference hits exceed u32 range".to_string())?;
    let mut previous = None;
    for _ in 0..hit_count {
        let hit = read_u32(input)?;
        if hit as usize >= reference_count {
            return Err("dictionary reference id is out of range".to_string());
        }
        if previous.is_some_and(|last| hit <= last) {
            return Err("dictionary reference hits are not strictly increasing".to_string());
        }
        previous = Some(hit);
        if hit_count > 1 {
            packed_hits.push(hit);
        }
    }
    if hit_count == 1 {
        Ok(ReferenceHits::One(previous.expect("nonempty hits")))
    } else {
        let len = u32::try_from(hit_count)
            .map_err(|_| "reference hit list exceeds u32 range".to_string())?;
        offset
            .checked_add(len)
            .ok_or_else(|| "packed reference hit range exceeds u32 range".to_string())?;
        Ok(ReferenceHits::Packed { offset, len })
    }
}

// 每次落盘的文件句柄留着复用,免得成千上万次 flush 都重新 open/close 一遍。
#[cfg(unix)]
fn raise_fd_limit() {
    unsafe {
        let mut limit: libc::rlimit = std::mem::zeroed();
        if libc::getrlimit(libc::RLIMIT_NOFILE, &mut limit) == 0 && limit.rlim_cur < limit.rlim_max
        {
            let mut raised = limit;
            raised.rlim_cur = limit.rlim_max;
            let _ = libc::setrlimit(libc::RLIMIT_NOFILE, &raised);
        }
    }
}

#[cfg(not(unix))]
fn raise_fd_limit() {}

struct OutputManager {
    paths: Vec<PathBuf>,
    buffers: Vec<Vec<u8>>,
    handles: Vec<Option<File>>,
    total_buffered: usize,
    file_budget: usize,
    paired_paths: bool,
    binary: bool,
    encode_scratch: Vec<u8>,
    recycled_buffers: Vec<Vec<u8>>,
}

impl OutputManager {
    fn new(
        output: &Path,
        subdir: &str,
        names: &[String],
        kind: FileKind,
        mode: u8,
    ) -> AppResult<Self> {
        let directory = output.join(subdir);
        fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
        let extension = if kind == FileKind::Fasta {
            ".fasta"
        } else {
            ".fq"
        };
        let mut paths = Vec::new();
        match mode {
            0 => {
                for name in names {
                    paths.push(directory.join(format!("{name}{extension}")));
                }
            }
            1 => {
                paths.push(directory.join("all_1.fq"));
                paths.push(directory.join("all_2.fq"));
            }
            3 => {}
            4 => {
                for name in names {
                    paths.push(directory.join(format!("{name}_1{extension}")));
                    paths.push(directory.join(format!("{name}_2{extension}")));
                }
            }
            5 => {
                for name in names {
                    paths.push(directory.join(format!("{name}_1.gm2")));
                    paths.push(directory.join(format!("{name}_2.gm2")));
                }
            }
            _ => return Err("output manager does not support this mode".to_string()),
        }
        for path in &paths {
            if path.exists() {
                fs::remove_file(path).map_err(|e| e.to_string())?;
            }
        }
        let file_budget = if paths.is_empty() {
            MIN_FILE_BUDGET
        } else {
            (TOTAL_BUFFER_BUDGET / paths.len()).max(MIN_FILE_BUDGET)
        };
        let buffer_count = paths.len();
        Ok(Self {
            paths,
            buffers: vec![Vec::new(); buffer_count],
            handles: (0..buffer_count).map(|_| None).collect(),
            total_buffered: 0,
            file_budget,
            paired_paths: mode == 4 || mode == 5,
            binary: mode == 5,
            encode_scratch: Vec::with_capacity(4096),
            recycled_buffers: Vec::with_capacity(16),
        })
    }

    fn write_encoded_record(
        &mut self,
        reference: usize,
        encoded: &[u8],
        second_mate: bool,
    ) -> AppResult<()> {
        let key = if self.paired_paths {
            2 * reference + usize::from(second_mate)
        } else {
            reference
        };
        self.write_encoded_key(key, encoded)
    }

    fn write_key(&mut self, key: usize, record: &Record) -> AppResult<()> {
        if self.paths.is_empty() {
            return Ok(());
        }
        encode_record_into(record, self.binary, &mut self.encode_scratch)?;
        let encoded = std::mem::take(&mut self.encode_scratch);
        let result = self.write_encoded_key(key, &encoded);
        self.encode_scratch = encoded;
        result
    }

    fn write_encoded_key(&mut self, key: usize, encoded: &[u8]) -> AppResult<()> {
        if self.paths.is_empty() {
            return Ok(());
        }
        if self.buffers[key].capacity() == 0 {
            if let Some(buffer) = self.recycled_buffers.pop() {
                self.buffers[key] = buffer;
            }
        }
        self.total_buffered += encoded.len();
        self.buffers[key].extend_from_slice(encoded);
        if self.buffers[key].len() >= self.file_budget {
            self.flush_key(key)?;
        } else if self.total_buffered >= TOTAL_BUFFER_BUDGET {
            self.flush_to_low_watermark()?;
        }
        Ok(())
    }

    fn flush_key(&mut self, key: usize) -> AppResult<()> {
        if self.buffers[key].is_empty() {
            return Ok(());
        }
        if self.handles[key].is_none() {
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.paths[key])
                .map_err(|e| e.to_string())?;
            self.handles[key] = Some(file);
        }
        let out = self.handles[key].as_mut().expect("handle just populated");
        out.write_all(&self.buffers[key])
            .map_err(|e| e.to_string())?;
        let mut buffer = std::mem::take(&mut self.buffers[key]);
        self.total_buffered -= buffer.len();
        buffer.clear();
        // 仅回收中等大小的缓冲；异常大 locus 的峰值容量立即释放，避免常驻内存膨胀。
        const MAX_RECYCLED_BUFFER_CAPACITY: usize = 1024 * 1024;
        if buffer.capacity() <= MAX_RECYCLED_BUFFER_CAPACITY && self.recycled_buffers.len() < 16 {
            self.recycled_buffers.push(buffer);
        }
        Ok(())
    }

    // 越过 64 MiB 高水位时，一次挑出最大的活动缓冲，写到 32 MiB 以下。
    // 这避免均匀分散到大量 loci 时反复全表扫描、每次只释放一个小缓冲。
    fn flush_to_low_watermark(&mut self) -> AppResult<()> {
        let mut keys: Vec<usize> = self
            .buffers
            .iter()
            .enumerate()
            .filter_map(|(key, buffer)| (!buffer.is_empty()).then_some(key))
            .collect();
        keys.sort_unstable_by_key(|&key| std::cmp::Reverse(self.buffers[key].len()));
        for key in keys {
            if self.total_buffered <= BUFFER_LOW_WATERMARK {
                break;
            }
            self.flush_key(key)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> AppResult<()> {
        for key in 0..self.paths.len() {
            self.flush_key(key)?;
        }
        Ok(())
    }
}

// 文本输出保留原样，方便人眼复查。
fn encode_text_into(record: &Record, output: &mut Vec<u8>) {
    output.clear();
    output.extend_from_slice(&record.header);
    output.push(b'\n');
    output.extend_from_slice(&record.sequence);
    output.push(b'\n');
    if record.has_quality {
        output.extend_from_slice(&record.plus);
        output.push(b'\n');
        output.extend_from_slice(&record.quality);
        output.push(b'\n');
    }
}

fn encode_record_into(record: &Record, binary: bool, output: &mut Vec<u8>) -> AppResult<()> {
    if binary {
        encode_gm2_into(record, output)
    } else {
        encode_text_into(record, output);
        Ok(())
    }
}

fn gm2_base_value(base: u8) -> u8 {
    let upper = base.to_ascii_uppercase();
    if upper.is_ascii_uppercase() {
        upper - 64
    } else {
        b'N' - 64
    }
}

fn append_gm2_sequence_chunk(output: &mut Vec<u8>, last_chunk: &mut u8, chunk: u8) {
    let delta = chunk ^ *last_chunk;
    if output.len() > 6 && output.last().copied() == Some(delta) {
        if let Some(last) = output.last_mut() {
            *last |= 0x80;
        }
    } else {
        output.push(delta);
    }
    *last_chunk = chunk;
}

// GM2 把序列压紧，磁盘慢的时候少搬点儿字节就挺顶用。
fn encode_gm2_into(record: &Record, output: &mut Vec<u8>) -> AppResult<()> {
    output.clear();
    if record.sequence.is_empty() {
        return Ok(());
    }
    if record.sequence.len() > 0x7f_ffff {
        return Err("GM2 sequence is too long".to_string());
    }
    let capacity = record
        .sequence
        .len()
        .checked_mul(2)
        .and_then(|size| size.checked_add(6))
        .ok_or_else(|| "GM2 record size overflow".to_string())?;
    if output.capacity() < capacity {
        output.reserve(capacity - output.capacity());
    }
    output.resize(6, 0);

    let mut last_chunk = 0_u8;
    let mut last_value = gm2_base_value(record.sequence[0]);
    let mut duplicate_count = 0_u8;
    for &base in &record.sequence[1..] {
        let value = gm2_base_value(base);
        if value != last_value || duplicate_count == 3 {
            append_gm2_sequence_chunk(output, &mut last_chunk, (duplicate_count << 5) | last_value);
            last_value = value;
            duplicate_count = 0;
        } else {
            duplicate_count += 1;
        }
    }
    append_gm2_sequence_chunk(output, &mut last_chunk, (duplicate_count << 5) | last_value);

    let has_quality = record.has_quality && record.quality.len() >= record.sequence.len();
    if has_quality {
        let quality = &record.quality;
        let mut start = 0;
        while start < record.sequence.len() {
            let value = quality[start] & 0x7f;
            let mut length = 1_usize;
            while start + length < record.sequence.len()
                && quality[start + length] & 0x7f == value
                && length < 256
            {
                length += 1;
            }
            if length == 1 {
                output.push(value);
            } else {
                let duplicate_count = length - 1;
                output.push(0x80 | (duplicate_count as u8 & 0x7f));
                output.push(value | if duplicate_count & 0x80 != 0 { 0x80 } else { 0 });
            }
            start += length;
        }
    }
    let payload_length = output.len() - 6;
    if payload_length > 0xff_ffff {
        return Err("GM2 payload is too long".to_string());
    }
    output[..6].copy_from_slice(&[
        ((payload_length >> 16) & 0xff) as u8,
        ((payload_length >> 8) & 0xff) as u8,
        (payload_length & 0xff) as u8,
        (((record.sequence.len() >> 16) & 0x7f) as u8) | if has_quality { 0x80 } else { 0 },
        ((record.sequence.len() >> 8) & 0xff) as u8,
        (record.sequence.len() & 0xff) as u8,
    ]);
    Ok(())
}

#[cfg(test)]
fn encode_gm2(record: &Record) -> AppResult<Vec<u8>> {
    let mut output = Vec::new();
    encode_gm2_into(record, &mut output)?;
    Ok(output)
}

struct Logger {
    file: BufWriter<File>,
    last_flush: Instant,
}

impl Logger {
    fn new(output: &Path) -> AppResult<Self> {
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(output.join("log.txt"))
            .map_err(|e| e.to_string())?;
        Ok(Self {
            file: BufWriter::with_capacity(64 * 1024, file),
            last_flush: Instant::now(),
        })
    }
    fn log(&mut self, message: &str) {
        println!("{message}");
        let _ = writeln!(self.file, "{message}");
        if self.last_flush.elapsed().as_millis() >= 500 {
            let _ = self.file.flush();
            self.last_flush = Instant::now();
        }
    }
}
impl Drop for Logger {
    fn drop(&mut self) {
        let _ = self.file.flush();
    }
}

// 每个 locus 最后留个计数单子，方便看过滤到底捞着多少 reads。
// 仅扫描模式是兼容/诊断功能：保持流式单线程，不保存批次也不创建输出管理器。
fn scan_mode3_single(
    args: &Args,
    index: &KmerIndex,
    fallbacks: &[KmerIndex],
    siblings: &[Option<u32>],
    kind: FileKind,
    logger: &mut Logger,
) -> AppResult<Vec<u64>> {
    let reference_count = index.reference_names.len();
    let mut counts = vec![0_u64; reference_count];
    let mut collector = HitCollector::new(reference_count);
    collector.minimum_span = args.minimum_hit_span;
    let mut direct_counts = vec![0_u64; reference_count];
    let max_reads = args.max_read_blocks.saturating_mul(MEBIBYTE_READS);
    for file_number in 0..args.q1.len() {
        let mut reader1 = SequenceReader::open(&args.q1[file_number], kind, false)?;
        let mut reader2 = if args.q2.is_empty() {
            None
        } else {
            Some(SequenceReader::open(&args.q2[file_number], kind, false)?)
        };
        let mut record1 = Record::default();
        let mut record2 = Record::default();
        let mut read_count = 0_u64;
        let mut stopped_at_limit = false;
        while reader1.next_record_into(&mut record1)? {
            let paired = match reader2.as_mut() {
                Some(reader) => {
                    if !reader.next_record_into(&mut record2)? {
                        return Err(
                            "paired input files contain different numbers of records".to_string()
                        );
                    }
                    true
                }
                None => false,
            };
            collector.begin();
            collect_fragment_hits(
                index,
                fallbacks,
                &record1.sequence,
                paired.then_some(record2.sequence.as_slice()),
                args.step,
                &mut collector,
            );
            if args.link_rad_arms {
                collector.link_siblings(siblings, &mut direct_counts, args.link_rad_max_fragments);
            }
            let increment = if paired { 2 } else { 1 };
            for &reference in &collector.hits {
                counts[reference as usize] += increment;
            }
            read_count += 1;
            if read_count.is_multiple_of(MEBIBYTE_READS) {
                logger.log(&format!(
                    "Handled {} Mi reads from {}.",
                    read_count / MEBIBYTE_READS,
                    args.q1[file_number].display()
                ));
            }
            if args.max_read_blocks > 0 && read_count >= max_reads {
                stopped_at_limit = true;
                break;
            }
        }
        if !stopped_at_limit {
            if let Some(reader) = reader2.as_mut() {
                if reader.next_record_into(&mut record2)? {
                    return Err(
                        "paired input files contain different numbers of records".to_string()
                    );
                }
            }
        }
    }
    Ok(counts)
}

fn write_read_counts(output: &Path, names: &[String], counts: &[u64]) -> AppResult<()> {
    let mut out =
        File::create(output.join("ref_reads_count_dict.txt")).map_err(|e| e.to_string())?;
    for (name, count) in names.iter().zip(counts) {
        if *count > 0 {
            writeln!(out, "{name},{count}").map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn run(args: Args) -> AppResult<()> {
    // 索引只建一次，随后按输入 read 逐条分到命中的 locus。
    fs::create_dir_all(&args.output).map_err(|e| e.to_string())?;
    let mut logger = Logger::new(&args.output)?;
    logger.log("Getting information from references...");
    let started = Instant::now();
    let cached = args
        .dictionary
        .as_ref()
        .filter(|path| path.exists())
        .map(|path| {
            probe_dictionary(path, args.kmer)?;
            let reference_hash = reference_content_hash(&args.reference)?;
            load_dictionary(path, args.kmer, &reference_hash)
        });
    let index = match cached {
        Some(Ok(index)) => {
            logger.log(&format!(
                "Loaded k-mer dictionary with {} entries.",
                index.len()
            ));
            index
        }
        Some(Err(reason)) => {
            logger.log(&format!("Ignoring reusable dictionary: {reason}"));
            let (index, hash) = build_index(&args.reference, args.kmer, true)?;
            logger.log(&format!(
                "Built k-mer dictionary with {} entries.",
                index.len()
            ));
            if let (Some(path), Some(hash)) = (&args.dictionary, hash.as_ref()) {
                write_dictionary(&index, hash, path)?;
            }
            index
        }
        None => {
            let need_hash = args.dictionary.is_some();
            let (index, hash) = build_index(&args.reference, args.kmer, need_hash)?;
            logger.log(&format!(
                "Built k-mer dictionary with {} entries.",
                index.len()
            ));
            if let (Some(path), Some(hash)) = (&args.dictionary, hash.as_ref()) {
                write_dictionary(&index, hash, path)?;
            }
            index
        }
    };
    let mut fallback_indices = Vec::with_capacity(args.fallback_kmers.len());
    for &fallback_k in &args.fallback_kmers {
        let (fallback, _) = build_index(&args.reference, fallback_k, false)?;
        logger.log(&format!(
            "Built fallback {fallback_k}-mer dictionary with {} entries.",
            fallback.len()
        ));
        fallback_indices.push(fallback);
    }
    let siblings = if args.link_rad_arms {
        rad_arm_siblings(&index.reference_names)?
    } else {
        vec![None; index.reference_names.len()]
    };
    logger.log(&format!(
        "Dictionary stage took {:.3} seconds.",
        started.elapsed().as_secs_f64()
    ));
    if args.mode == 2 {
        return Ok(());
    }
    if args.use_composition_pattern {
        logger
            .log("Note: -lb is accepted for compatibility; exact rolling k-mer scanning is used.");
    }

    let kind = detect_kind(&args.q1[0])?;
    for path in args.q1.iter().chain(args.q2.iter()) {
        if detect_kind(path)? != kind {
            return Err("all read files must use the same FASTA/FASTQ format".to_string());
        }
    }
    if args
        .q1
        .iter()
        .chain(args.q2.iter())
        .any(|path| is_gzip(path))
    {
        logger.log(&format!("Gzip backend: {}.", zlib_backend().name));
    }
    if args.mode == 3 {
        let filter_started = Instant::now();
        let counts = scan_mode3_single(
            &args,
            &index,
            &fallback_indices,
            &siblings,
            kind,
            &mut logger,
        )?;
        write_read_counts(&args.output, &index.reference_names, &counts)?;
        logger.log(&format!(
            "Filtering took {:.3} seconds.",
            filter_started.elapsed().as_secs_f64()
        ));
        return Ok(());
    }
    let keep_text_lines = matches!(args.mode, 0 | 1 | 4);
    let filter_started = Instant::now();
    let mut output = OutputManager::new(
        &args.output,
        &args.out_subdir,
        &index.reference_names,
        kind,
        args.mode,
    )?;
    let mut counts = vec![0_u64; index.reference_names.len()];
    let mut direct_counts = vec![0_u64; index.reference_names.len()];
    let mut collector = HitCollector::new(index.reference_names.len());
    collector.minimum_span = args.minimum_hit_span;
    let mut encoded1 = Vec::with_capacity(4096);
    let mut encoded2 = Vec::with_capacity(4096);
    if args.get_reverse {
        logger.log(
            "Note: -gr is retained for compatibility; canonical k-mers already match both strands.",
        );
    }
    for file_number in 0..args.q1.len() {
        let mut reader1 = SequenceReader::open(&args.q1[file_number], kind, keep_text_lines)?;
        let mut reader2 = if args.q2.is_empty() {
            None
        } else {
            Some(SequenceReader::open(
                &args.q2[file_number],
                kind,
                keep_text_lines,
            )?)
        };
        let mut read_count = 0_u64;
        let max_reads = args.max_read_blocks.saturating_mul(MEBIBYTE_READS);
        let mut stopped_at_limit = false;
        let mut record1 = Record::default();
        let mut record2 = Record::default();
        while reader1.next_record_into(&mut record1)? {
            let has_record2 = match reader2.as_mut() {
                Some(reader) => {
                    if !reader.next_record_into(&mut record2)? {
                        return Err(
                            "paired input files contain different numbers of records".to_string()
                        );
                    }
                    true
                }
                None => false,
            };
            collector.begin();
            collect_fragment_hits(
                &index,
                &fallback_indices,
                &record1.sequence,
                has_record2.then_some(record2.sequence.as_slice()),
                args.step,
                &mut collector,
            );
            if args.link_rad_arms {
                collector.link_siblings(&siblings, &mut direct_counts, args.link_rad_max_fragments);
            }
            if !collector.hits.is_empty() {
                if args.mode == 1 {
                    output.write_key(0, &record1)?;
                    if has_record2 {
                        output.write_key(1, &record2)?;
                    }
                }
                let write_per_reference = matches!(args.mode, 0 | 4 | 5);
                if write_per_reference {
                    encode_record_into(&record1, args.mode == 5, &mut encoded1)?;
                    if has_record2 {
                        encode_record_into(&record2, args.mode == 5, &mut encoded2)?;
                    }
                }
                for &reference in &collector.hits {
                    let reference = reference as usize;
                    counts[reference] += if has_record2 { 2 } else { 1 };
                    if write_per_reference {
                        output.write_encoded_record(reference, &encoded1, false)?;
                    }
                    if write_per_reference && has_record2 {
                        output.write_encoded_record(reference, &encoded2, true)?;
                    }
                }
            }
            read_count += 1;
            if read_count.is_multiple_of(MEBIBYTE_READS) {
                logger.log(&format!(
                    "Handled {} Mi reads from {}.",
                    read_count / MEBIBYTE_READS,
                    args.q1[file_number].display()
                ));
            }
            if args.max_read_blocks > 0 && read_count >= max_reads {
                stopped_at_limit = true;
                break;
            }
        }
        if !stopped_at_limit {
            if let Some(reader) = reader2.as_mut() {
                let mut probe = Record::default();
                if reader.next_record_into(&mut probe)? {
                    return Err(
                        "paired input files contain different numbers of records".to_string()
                    );
                }
            }
        }
    }

    output.flush()?;
    write_read_counts(&args.output, &index.reference_names, &counts)?;
    logger.log(&format!(
        "Filtering took {:.3} seconds.",
        filter_started.elapsed().as_secs_f64()
    ));
    Ok(())
}

fn main() {
    // 输出文件句柄常驻复用,先把 fd 上限尽量抬到硬上限,避免 locus 数一多就撞到 "too many open files"。
    raise_fd_limit();
    let argv: Vec<String> = env::args().skip(1).collect();
    let args = match parse_args(argv) {
        Ok(args) => args,
        Err(error) => {
            eprintln!("Invalid argument: {error}");
            process::exit(2);
        }
    };
    if let Err(error) = run(args) {
        eprintln!("Error: {error}");
        process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn normalized_kmer(sequence: &[u8]) -> Option<Vec<u8>> {
        sequence
            .iter()
            .map(|&base| base_code(base).map(|code| b"ACGT"[code as usize]))
            .collect()
    }

    fn reverse_complement(sequence: &[u8]) -> Option<Vec<u8>> {
        sequence
            .iter()
            .rev()
            .map(|&base| base_code(base).map(|code| b"ACGT"[(3 - code) as usize]))
            .collect()
    }

    fn oracle_hits(
        references: &[Vec<u8>],
        read: &[u8],
        k: usize,
        step: usize,
        reverse_indexed: bool,
    ) -> Vec<u32> {
        if read.len() < k {
            return Vec::new();
        }
        let reference_kmers: Vec<HashSet<Vec<u8>>> = references
            .iter()
            .map(|reference| {
                let mut kmers = HashSet::new();
                for window in reference.windows(k) {
                    if let Some(forward) = normalized_kmer(window) {
                        if reverse_indexed {
                            kmers.insert(reverse_complement(&forward).unwrap());
                        }
                        kmers.insert(forward);
                    }
                }
                kmers
            })
            .collect();
        let tail = read.len() - k;
        let mut starts: Vec<usize> = (0..=tail).step_by(step).collect();
        if starts.last().copied() != Some(tail) {
            starts.push(tail);
        }
        let mut hits = Vec::new();
        for (reference, kmers) in reference_kmers.iter().enumerate() {
            let matched = starts.iter().any(|&start| {
                let window = &read[start..start + k];
                normalized_kmer(window).is_some_and(|forward| {
                    kmers.contains(&forward)
                        || (!reverse_indexed
                            && reverse_complement(&forward)
                                .is_some_and(|reverse| kmers.contains(&reverse)))
                })
            });
            if matched {
                hits.push(reference as u32);
            }
        }
        hits
    }

    fn pseudo_random_sequence(seed: &mut u64, length: usize) -> Vec<u8> {
        let alphabet = b"ACGTACGTACGTACGTACGTACGTACGTACGTUN";
        (0..length)
            .map(|_| {
                *seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                alphabet[(*seed as usize) % alphabet.len()]
            })
            .collect()
    }

    #[test]
    fn sensitive_recruitment_requires_a_span_or_independent_mates_per_family() {
        let sequence = b"AGCTTAGCGACCTGACTGACGTATCAGTA";
        let unrelated = b"CGATTGCAATGCCATAGTCATACGGTCAA";
        let mut index = KmerIndex::new(21, vec!["a".into(), "b".into()]);
        index.add_reference_sequence(sequence, 0);
        index.add_reference_sequence(unrelated, 1);
        index.finalize_hits().unwrap();
        let mut hits = HitCollector::new(2);
        hits.minimum_span = 28;
        hits.begin();
        collect_fragment_hits(&index, &[], &sequence[..21], None, 4, &mut hits);
        assert!(hits.hits.is_empty());
        hits.begin();
        collect_fragment_hits(&index, &[], sequence, None, 4, &mut hits);
        assert_eq!(hits.hits, vec![0]);
        hits.begin();
        collect_fragment_hits(
            &index,
            &[],
            &sequence[..21],
            Some(&sequence[..21]),
            4,
            &mut hits,
        );
        assert_eq!(hits.hits, vec![0]);
        hits.begin();
        collect_fragment_hits(
            &index,
            &[],
            &sequence[..21],
            Some(&unrelated[..21]),
            4,
            &mut hits,
        );
        assert!(hits.hits.is_empty());
    }

    #[test]
    fn reference_hits_layout_is_already_compact() {
        assert_eq!(std::mem::size_of::<ReferenceHits>(), 12);
    }

    #[test]
    fn ambiguous_bases_split_reference_kmers() {
        let mut index = KmerIndex::new(16, vec!["locus".to_string()]);
        index.add_reference_sequence(b"AAAAAAAAAAAAAAAANCCCCCCCCCCCCCCCC", 0);
        assert_eq!(index.len(), 2);
        let KmerStore::Short(map) = &index.store else {
            unreachable!();
        };
        assert!(!map.contains_key(&0b01));
    }

    #[test]
    fn reverse_complement_lookup_matches_original_modes() {
        let mut index = KmerIndex::new(16, vec!["locus".to_string()]);
        index.add_reference_sequence(b"AAAACCCCGGGGTTTA", 0);
        let mut hits = HitCollector::new(1);
        hits.begin();
        index.collect_hits(b"TAAACCCCGGGGTTTT", 3, &mut hits);
        assert_eq!(hits.hits, vec![0]);
    }

    #[test]
    fn dictionary_roundtrip_preserves_medium_kmers_and_packed_hits() {
        let mut index = KmerIndex::new(33, vec!["a".to_string(), "b".to_string()]);
        index.add_reference_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT", 0);
        index.add_reference_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT", 1);
        index.finalize_hits().unwrap();
        let KmerStore::Medium(_) = &index.store else {
            panic!("33-mer must use u128 store")
        };
        let path = std::env::temp_dir().join(format!(
            "gm2-mainfilter-{}-{}.dict",
            process::id(),
            Instant::now().elapsed().as_nanos()
        ));
        let hash = [7_u8; 32];
        write_dictionary(&index, &hash, &path).unwrap();
        let loaded = load_dictionary(&path, 33, &hash).unwrap();
        let _ = fs::remove_file(path);
        let KmerStore::Medium(map) = &loaded.store else {
            panic!("loaded 33-mer must use u128 store")
        };
        assert_eq!(map.len(), index.len());
        let hits = map.values().next().unwrap();
        assert_eq!(loaded.hit_slice(hits), &[0, 1]);
    }

    #[test]
    fn gm2_header_uses_all_six_bytes() {
        let record = Record {
            sequence: vec![b'A'; 300],
            quality: vec![b'I'; 300],
            has_quality: true,
            ..Record::default()
        };
        let encoded = encode_gm2(&record).unwrap();
        let payload =
            ((encoded[0] as usize) << 16) | ((encoded[1] as usize) << 8) | encoded[2] as usize;
        let sequence = (((encoded[3] & 0x7f) as usize) << 16)
            | ((encoded[4] as usize) << 8)
            | encoded[5] as usize;
        assert_eq!(encoded.len(), payload + 6);
        assert_eq!(sequence, 300);
        assert_ne!(encoded[3] & 0x80, 0);
    }

    #[test]
    fn rad_arm_links_are_bidirectional_and_locus_scoped() {
        let names = vec![
            "locus_a__R1".into(),
            "locus_b__R1".into(),
            "locus_a__R2".into(),
            "unpaired".into(),
        ];
        let siblings = rad_arm_siblings(&names).unwrap();
        assert_eq!(siblings, vec![Some(2), None, Some(0), None]);
        let mut hits = HitCollector::new(names.len());
        hits.begin();
        hits.mark(0);
        hits.link_siblings(&siblings, &mut vec![0; names.len()], 10_000);
        assert_eq!(hits.hits, vec![0, 2]);
    }

    #[test]
    fn rad_arm_linking_stops_at_fragment_cap() {
        let siblings = vec![Some(1), Some(0)];
        let mut direct_counts = vec![0; 2];
        let mut hits = HitCollector::new(2);

        hits.begin();
        hits.mark(0);
        hits.link_siblings(&siblings, &mut direct_counts, 1);
        assert_eq!(hits.hits, vec![0, 1]);

        hits.begin();
        hits.mark(0);
        hits.link_siblings(&siblings, &mut direct_counts, 1);
        assert_eq!(hits.hits, vec![0]);
    }

    #[test]
    fn fallback_index_is_queried_only_after_primary_misses() {
        let names = vec!["a".into(), "b".into()];
        let sequences = [
            b"ACGTCAGTGCATGACTCAGTACGA".as_slice(),
            b"TTGCAAGCTTAGGCTAACCGTTAA".as_slice(),
        ];
        let mut primary = KmerIndex::new(16, names.clone());
        let mut fallback = KmerIndex::new(8, names);
        for (reference, sequence) in sequences.iter().enumerate() {
            primary.add_reference_sequence(sequence, reference as u32);
            fallback.add_reference_sequence(sequence, reference as u32);
        }
        primary.finalize_hits().unwrap();
        fallback.finalize_hits().unwrap();
        let mut hits = HitCollector::new(2);
        hits.begin();
        collect_fragment_hits(
            &primary,
            std::slice::from_ref(&fallback),
            sequences[0],
            None,
            1,
            &mut hits,
        );
        assert_eq!(hits.hits, vec![0]);

        hits.begin();
        collect_fragment_hits(
            &primary,
            &[fallback],
            b"TTGCAAGCAAAAAAAAAAAAAAAA",
            None,
            1,
            &mut hits,
        );
        assert_eq!(hits.hits, vec![1]);
    }

    #[test]
    fn randomized_lookup_matches_naive_oracle_across_boundaries() {
        let mut seed = 0x4d59_5df4_d0f3_3173;
        for &k in &[16, 31, 32, 33, 48, 64, 65] {
            for &step in &[1, 3, k + 2] {
                for reverse_indexed in [false, true] {
                    for iteration in 0..8 {
                        let references: Vec<Vec<u8>> = (0..3)
                            .map(|_| pseudo_random_sequence(&mut seed, 96))
                            .collect();
                        let mut read = if iteration % 3 == 0 {
                            pseudo_random_sequence(&mut seed, 64)
                        } else {
                            references[iteration % references.len()][12..76].to_vec()
                        };
                        if iteration % 2 == 1 {
                            read = reverse_complement(&read).unwrap_or(read);
                        }
                        if iteration % 4 == 2 {
                            read[3] = b'N';
                        }

                        let mut index = KmerIndex::new(
                            k,
                            vec!["a".to_string(), "b".to_string(), "c".to_string()],
                        );
                        for (reference, sequence) in references.iter().enumerate() {
                            index.add_reference_sequence(sequence, reference as u32);
                        }
                        let mut collector = HitCollector::new(references.len());
                        collector.begin();
                        index.collect_hits(&read, step, &mut collector);
                        collector.hits.sort_unstable();

                        assert_eq!(
                            collector.hits,
                            oracle_hits(&references, &read, k, step, reverse_indexed),
                            "k={k}, step={step}, reverse_indexed={reverse_indexed}, iteration={iteration}"
                        );
                    }
                }
            }
        }
    }
}
