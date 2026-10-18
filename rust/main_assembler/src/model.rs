use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum AssemblyMode {
    Reference,
    Uce,
    Its2,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum PathStrategy {
    Search,
    Backbone,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GraphFormat {
    None,
    Gfa,
    Dot,
    Both,
}

#[derive(Clone, Debug)]
pub struct Args {
    pub reference: PathBuf,
    pub output: PathBuf,
    pub kmer_size: usize,
    pub kmer_min: usize,
    pub kmer_max: usize,
    pub error_limit: u32,
    pub iteration: usize,
    pub min_coverage: f64,
    pub soft_boundary: i64,
    pub threads: usize,
    pub assembly_mode: AssemblyMode,
    pub side_candidates: usize,
    pub path_strategy: PathStrategy,
    pub backbone_lookahead: usize,
    pub reverse_reuse_reference_scale: f64,
    pub max_contig_length: usize,
    pub min_read_density: f64,
    pub density_check_min_length: usize,
    pub max_depth_cv: f64,
    pub read_chunk_size: usize,
    pub kmer_count_threads: usize,
    pub graph_format: GraphFormat,
    pub max_depth_ratio: f64,
    pub reference_cache_dir: Option<PathBuf>,
    pub profile: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct RefKmer {
    pub depth: u32,
    pub position: i32,
    pub is_reverse: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct KmerInfo {
    pub depth: i64,
    pub position: i32,
    pub is_reverse: bool,
    pub reference_weight: i64,
    pub fragment_support: u32,
    pub pe_flank_rescue: bool,
}

#[derive(Clone, Copy, Debug)]
pub struct Node {
    pub kmer: u128,
    pub position: i32,
    pub weight: i64,
    pub pe_support: u32,
}

#[derive(Clone, Debug, Default)]
pub struct PathContig {
    pub weights: Vec<i64>,
    pub bases: Vec<u8>,
}

#[derive(Clone, Debug, Default)]
pub struct SideContig {
    pub sequence: Vec<u8>,
    pub weight: i64,
    pub read_count: u64,
}

#[derive(Clone, Debug, Default)]
pub struct ReadSupport {
    pub total_read_count: u64,
    pub unique_read_count: u64,
    pub multi_mapping_read_count: u64,
    pub supported_extent: usize,
    pub supported_bases: usize,
    pub breadth: f64,
    pub max_gap: usize,
    pub flank_balance: f64,
    pub left_coord: usize,
    pub right_coord: usize,
}

#[derive(Clone, Debug, Default)]
pub struct ContigRecord {
    pub label: String,
    pub equivalence_members: String,
    pub sequence: Vec<u8>,
    pub its2_full_reference_span: bool,
    pub seed_count: usize,
    pub position: i32,
    pub weight: i64,
    pub read_count: u64,
    pub supported_span: usize,
    pub flank_balance: f64,
    pub read_density: f64,
    pub support_fraction: f64,
    pub kmer_median_depth: f64,
    pub kmer_depth_cv: f64,
    pub kmer_max_depth_ratio: f64,
    pub unique_read_count: u64,
    pub multi_mapping_read_count: u64,
    pub supported_bases: usize,
    pub support_breadth: f64,
    pub max_support_gap: usize,
    pub fragment_support: u64,
    pub paired_fragment_support: u64,
    pub diagnostic_fragment_support: u64,
    pub quality_rescued_fragment_support: u64,
    pub phase_components: Option<usize>,
    pub min_variant_site_support: Option<usize>,
    pub em_fragment_support: f64,
    pub em_abundance: f64,
    pub accepted: bool,
    pub rejection_reason: String,
}

#[derive(Clone, Debug)]
pub struct LocusTask {
    pub key: String,
    pub reference_path: PathBuf,
    pub reference_count: usize,
    pub ordinal: usize,
    pub total: usize,
}

#[derive(Clone, Debug, Default)]
pub struct LocusResult {
    pub key: String,
    pub status: String,
    pub value: u64,
    pub accepted: bool,
    pub rejection_reason: String,
    pub selected_contig_length: usize,
    pub read_supported_span: usize,
    pub slice_supported_bases: usize,
    pub slice_support_breadth: f64,
    pub max_slice_support_gap: usize,
    pub read_count: u64,
    pub unique_read_count: u64,
    pub multi_mapping_read_count: u64,
    pub read_density: f64,
    pub unique_read_density: f64,
    pub support_fraction: f64,
    pub flank_balance: f64,
    pub kmer_median_depth: f64,
    pub kmer_depth_cv: f64,
    pub kmer_max_depth_ratio: f64,
    pub candidate_count: usize,
    pub low_quality: bool,
    pub skipped: bool,
}

impl LocusResult {
    pub fn failure(key: &str, status: &str) -> Self {
        Self {
            key: key.to_string(),
            status: status.to_string(),
            ..Self::default()
        }
    }
}
