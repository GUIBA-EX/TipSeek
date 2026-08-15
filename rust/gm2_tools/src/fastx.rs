//! Buffered byte-level FASTA/FASTQ readers shared by TipSeek tools.
//! Text is kept as bytes: sequence tools do not need UTF-8 validation or a
//! temporary `String` for every FASTQ line.
use std::ffi::{c_char, c_int, c_uint, c_void, CString};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::sync::OnceLock;

pub const READ_BUFFER_SIZE: usize = 1024 * 1024;

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

#[cfg(unix)]
/// Resolves a dynamic zlib-ng symbol to a caller-supplied function-pointer type.
///
/// # Safety
///
/// `handle` must be a live handle returned by `dlopen`, and `F` must exactly
/// match the ABI and signature exported for `symbol`. Callers keep the library
/// loaded for at least as long as they use the returned pointer.
unsafe fn dlsym_typed<F: Copy>(handle: *mut c_void, symbol: &str) -> Option<F> {
    let symbol = CString::new(symbol).ok()?;
    let pointer = libc::dlsym(handle, symbol.as_ptr());
    (!pointer.is_null()).then(|| std::mem::transmute_copy(&pointer))
}

#[cfg(unix)]
fn detect_zlib_ng() -> Option<ZlibBackend> {
    for library in ["libz-ng.so.2", "libz-ng.so.1", "libz-ng.so"] {
        let library = CString::new(library).expect("static string contains no NUL");
        let handle = unsafe { libc::dlopen(library.as_ptr(), libc::RTLD_NOW | libc::RTLD_LOCAL) };
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
            // Keep the dynamic library resident while its function pointers are used.
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

static ZLIB_BACKEND: OnceLock<ZlibBackend> = OnceLock::new();

fn zlib_backend() -> ZlibBackend {
    *ZLIB_BACKEND.get_or_init(|| {
        build_detected_zlib_ng_backend()
            .or_else(detect_zlib_ng)
            .unwrap_or_else(stock_zlib_backend)
    })
}

/// Gzip reader with a 1 MiB compressed-input buffer.  It uses zlib-ng when a
/// compatible system library is available, otherwise the platform zlib ABI.
struct GzipReader {
    handle: *mut c_void,
    backend: ZlibBackend,
}

impl GzipReader {
    fn open(path: &Path) -> io::Result<Self> {
        let backend = zlib_backend();
        let path = CString::new(path.as_os_str().as_encoded_bytes())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains a NUL byte"))?;
        let mode = CString::new("rb").expect("static string contains no NUL");
        let handle = unsafe { (backend.open)(path.as_ptr(), mode.as_ptr()) };
        if handle.is_null() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "cannot open gzip file",
            ));
        }
        unsafe { (backend.buffer)(handle, READ_BUFFER_SIZE as c_uint) };
        Ok(Self { handle, backend })
    }
}

impl Read for GzipReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        let result = unsafe {
            (self.backend.read)(
                self.handle,
                buffer.as_mut_ptr().cast(),
                buffer.len().min(c_int::MAX as usize) as u32,
            )
        };
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

/// Backend selected for gzip files in this process, for profiling/logging.
pub fn gzip_backend_name() -> &'static str {
    zlib_backend().name
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FastxFormat {
    Fasta,
    Fastq,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FastxRecord {
    /// Includes the leading `>` or `@` exactly as it appeared in the input.
    pub header: Vec<u8>,
    pub sequence: Vec<u8>,
    /// Empty for FASTA records; includes the leading `+` for FASTQ records.
    pub plus: Vec<u8>,
    /// Empty for FASTA records.
    pub quality: Vec<u8>,
}

pub struct FastxReader {
    input: BufReader<Box<dyn Read>>,
    format: FastxFormat,
    pending_header: Option<Vec<u8>>,
    scratch: Vec<u8>,
    finished: bool,
}

fn is_gzip(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
}

pub fn open_input(path: &Path) -> io::Result<BufReader<Box<dyn Read>>> {
    let input: Box<dyn Read> = if is_gzip(path) {
        Box::new(GzipReader::open(path)?)
    } else {
        Box::new(File::open(path)?)
    };
    Ok(BufReader::with_capacity(READ_BUFFER_SIZE, input))
}

impl FastxReader {
    pub fn open(path: &Path, format: FastxFormat) -> io::Result<Self> {
        Ok(Self {
            input: open_input(path)?,
            format,
            pending_header: None,
            scratch: Vec::with_capacity(512),
            finished: false,
        })
    }

    fn read_line(&mut self) -> io::Result<Option<Vec<u8>>> {
        self.scratch.clear();
        if self.input.read_until(b'\n', &mut self.scratch)? == 0 {
            return Ok(None);
        }
        while matches!(self.scratch.last(), Some(b'\n' | b'\r')) {
            self.scratch.pop();
        }
        Ok(Some(std::mem::take(&mut self.scratch)))
    }

    pub fn next_record(&mut self) -> io::Result<Option<FastxRecord>> {
        match self.format {
            FastxFormat::Fasta => self.next_fasta(),
            FastxFormat::Fastq => self.next_fastq(),
        }
    }

    fn next_fastq(&mut self) -> io::Result<Option<FastxRecord>> {
        let header = loop {
            let Some(line) = self.read_line()? else {
                return Ok(None);
            };
            if !line.is_empty() {
                break line;
            }
        };
        if header.first() != Some(&b'@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed FASTQ record",
            ));
        }
        let sequence = self.read_line()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ sequence")
        })?;
        let plus = self.read_line()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ plus line")
        })?;
        if plus.first() != Some(&b'+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed FASTQ plus line",
            ));
        }
        let quality = self.read_line()?.ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ quality")
        })?;
        if quality.len() != sequence.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "FASTQ sequence and quality lengths differ",
            ));
        }
        Ok(Some(FastxRecord {
            header,
            sequence,
            plus,
            quality,
        }))
    }

    fn next_fasta(&mut self) -> io::Result<Option<FastxRecord>> {
        if self.finished {
            return Ok(None);
        }
        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                let Some(line) = self.read_line()? else {
                    self.finished = true;
                    return Ok(None);
                };
                if line.first() == Some(&b'>') {
                    break line;
                }
                if !line.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "FASTA sequence encountered before a header",
                    ));
                }
            }
        };
        let mut sequence = Vec::new();
        loop {
            let Some(line) = self.read_line()? else {
                self.finished = true;
                break;
            };
            if line.first() == Some(&b'>') {
                self.pending_header = Some(line);
                break;
            }
            sequence.extend_from_slice(&line);
        }
        Ok(Some(FastxRecord {
            header,
            sequence,
            plus: Vec::new(),
            quality: Vec::new(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn parses_fastq_bytes_and_crlf() {
        let path = std::env::temp_dir().join(format!("gm2-fastx-{}", std::process::id()));
        let mut file = File::create(&path).unwrap();
        file.write_all(b"@read/1\r\nAcGT\r\n+\r\n!!!!\r\n").unwrap();
        drop(file);
        let mut reader = FastxReader::open(&path, FastxFormat::Fastq).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.header, b"@read/1");
        assert_eq!(record.sequence, b"AcGT");
        assert_eq!(record.quality, b"!!!!");
        assert!(reader.next_record().unwrap().is_none());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn parses_gzip_fastq() {
        use flate2::{write::GzEncoder, Compression};
        let path =
            std::env::temp_dir().join(format!("gm2-fastx-gzip-{}.fq.gz", std::process::id()));
        let file = File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(b"@read\nACGT\n+\n!!!!\n").unwrap();
        encoder.finish().unwrap();
        let mut reader = FastxReader::open(&path, FastxFormat::Fastq).unwrap();
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"ACGT");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn parses_multiline_fasta() {
        let path = std::env::temp_dir().join(format!("gm2-fastx-fasta-{}", std::process::id()));
        let mut file = File::create(&path).unwrap();
        file.write_all(b">one\nAC\nGT\n>two\nNN\n").unwrap();
        drop(file);
        let mut reader = FastxReader::open(&path, FastxFormat::Fasta).unwrap();
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"ACGT");
        assert_eq!(reader.next_record().unwrap().unwrap().sequence, b"NN");
        assert!(reader.next_record().unwrap().is_none());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn rejects_truncated_and_malformed_fastx_without_panicking() {
        let root = std::env::temp_dir().join(format!(
            "gm2-fastx-invalid-{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        for (name, format, contents, expected) in [
            (
                "truncated_sequence.fq",
                FastxFormat::Fastq,
                b"@read\n".as_slice(),
                "truncated FASTQ sequence",
            ),
            (
                "bad_plus.fq",
                FastxFormat::Fastq,
                b"@read\nAC\n-\n!!\n".as_slice(),
                "malformed FASTQ plus line",
            ),
            (
                "short_quality.fq",
                FastxFormat::Fastq,
                b"@read\nAC\n+\n!\n".as_slice(),
                "FASTQ sequence and quality lengths differ",
            ),
            (
                "missing_header.fa",
                FastxFormat::Fasta,
                b"ACGT\n".as_slice(),
                "FASTA sequence encountered before a header",
            ),
        ] {
            let path = root.join(name);
            std::fs::write(&path, contents).unwrap();
            let mut reader = FastxReader::open(&path, format).unwrap();
            let error = reader.next_record().unwrap_err();
            assert!(error.to_string().contains(expected), "{error}");
        }
        std::fs::remove_dir_all(root).unwrap();
    }
}
