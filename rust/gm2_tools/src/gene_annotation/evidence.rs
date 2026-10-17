//! Engine-neutral evidence required by the common structural QC pipeline.
use super::{parse_protein_cigar, Interval, RawGeneModel};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentEngine {
    Miniprot,
    Exonerate,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentEvidence {
    pub engine: AlignmentEngine,
    pub query_len: usize,
    pub query_interval: Interval,
    pub alignment_score: Option<i64>,
    pub frameshifts: usize,
    pub inframe_stops: usize,
    pub large_query_gaps: usize,
    pub query_gap_aa: usize,
    /// Offsets into the spliced CDS for each aligned query residue; query-only
    /// gaps are None. The length is checked against the extracted CDS.
    pub coding_map: Option<(Vec<Option<usize>>, usize)>,
}

impl AlignmentEvidence {
    pub fn from_miniprot(raw: &RawGeneModel) -> Result<Option<Self>, String> {
        let Some(paf) = &raw.paf else { return Ok(None) };
        let mut frameshifts = paf.frameshifts;
        let mut cigar_frameshifts = 0;
        let mut large_query_gaps = 0;
        let mut query_gap_aa: usize = 0;
        let mut offsets = Vec::new();
        let mut coding_nt = 0;
        if let Some(cigar) = &paf.cigar {
            for operation in parse_protein_cigar(cigar)? {
                match operation.op {
                    'M' => {
                        offsets.extend((0..operation.length).map(|i| Some(coding_nt + 3 * i)));
                        coding_nt += 3 * operation.length;
                    }
                    'I' => offsets.extend(std::iter::repeat_n(None, operation.length)),
                    'D' => coding_nt += 3 * operation.length,
                    'U' | 'V' => {
                        offsets.push(Some(coding_nt));
                        coding_nt += 3;
                    }
                    _ => {}
                }
                if matches!(operation.op, 'F' | 'G') {
                    cigar_frameshifts += 1;
                }
                if operation.op == 'I' {
                    query_gap_aa = query_gap_aa
                        .checked_add(operation.length)
                        .ok_or("query_gap_overflow")?;
                }
                if operation.op == 'I' && operation.length >= 10 {
                    large_query_gaps += 1;
                }
            }
        }
        if query_gap_aa > paf.query_interval.len() {
            return Err("query_gap_exceeds_span".into());
        }
        frameshifts = frameshifts.max(cigar_frameshifts);
        Ok(Some(Self {
            engine: AlignmentEngine::Miniprot,
            query_len: paf.query_len,
            query_interval: paf.query_interval,
            alignment_score: paf.alignment_score,
            frameshifts,
            inframe_stops: paf.inframe_stops,
            large_query_gaps,
            query_gap_aa,
            coding_map: (cigar_frameshifts == 0 && offsets.len() == paf.query_interval.len())
                .then_some((offsets, coding_nt)),
        }))
    }
}
