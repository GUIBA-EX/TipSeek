//! Conservative residue coordinates between orthologous reference proteins.
//! A position is usable only if all minimum-cost alignments agree. Insertions,
//! repeat-shift alternatives and weak local anchors are left unmapped.
const INF: u32 = u32::MAX / 4;
const MAX_CELLS: usize = 4_000_000;

struct Band {
    start: Vec<usize>,
    end: Vec<usize>,
    offset: Vec<usize>,
    values: Vec<u32>,
}
impl Band {
    fn new(n: usize, m: usize, width: usize) -> Option<Self> {
        let mut start = Vec::new();
        let mut end = Vec::new();
        let mut offset = Vec::new();
        let mut cells = 0;
        for i in 0..=n {
            let a = i.saturating_sub(width).min(m + 1);
            let b = (i + width + 1).min(m + 1);
            start.push(a);
            end.push(b);
            offset.push(cells);
            cells += b.saturating_sub(a);
            if cells > MAX_CELLS {
                return None;
            }
        }
        Some(Self {
            start,
            end,
            offset,
            values: vec![INF; cells],
        })
    }
    fn get(&self, i: usize, j: usize) -> u32 {
        if i >= self.start.len() || j < self.start[i] || j >= self.end[i] {
            INF
        } else {
            self.values[self.offset[i] + j - self.start[i]]
        }
    }
    fn set(&mut self, i: usize, j: usize, value: u32) {
        self.values[self.offset[i] + j - self.start[i]] = value;
    }
}

fn align(a: &[u8], b: &[u8], width: usize) -> Option<(Band, Band, u32)> {
    let (n, m) = (a.len(), b.len());
    let mut forward = Band::new(n, m, width)?;
    for i in 0..=n {
        for j in forward.start[i]..forward.end[i] {
            let mut cost = if i == 0 && j == 0 { 0 } else { INF };
            if i > 0 {
                cost = cost.min(forward.get(i - 1, j) + 2);
            }
            if j > 0 {
                cost = cost.min(forward.get(i, j - 1) + 2);
            }
            if i > 0 && j > 0 {
                cost = cost.min(forward.get(i - 1, j - 1) + u32::from(a[i - 1] != b[j - 1]));
            }
            forward.set(i, j, cost);
        }
    }
    let cost = forward.get(n, m);
    if cost >= INF {
        return None;
    }
    let mut backward = Band::new(n, m, width)?;
    for i in (0..=n).rev() {
        for j in (backward.start[i]..backward.end[i]).rev() {
            let mut value = if i == n && j == m { 0 } else { INF };
            if i < n {
                value = value.min(backward.get(i + 1, j) + 2);
            }
            if j < m {
                value = value.min(backward.get(i, j + 1) + 2);
            }
            if i < n && j < m {
                value = value.min(backward.get(i + 1, j + 1) + u32::from(a[i] != b[j]));
            }
            backward.set(i, j, value);
        }
    }
    Some((forward, backward, cost))
}

#[derive(Debug)]
pub struct ResidueMap {
    pub positions: Vec<Option<usize>>,
    pub status: &'static str,
    pub cost: Option<u32>,
}

pub fn map(source: &str, backbone: &str) -> ResidueMap {
    let (a, b) = (source.as_bytes(), backbone.as_bytes());
    let rejected = |status| ResidueMap {
        positions: vec![None; a.len()],
        status,
        cost: None,
    };
    if a == b {
        return ResidueMap {
            positions: (0..a.len()).map(Some).collect(),
            status: "identical_reference",
            cost: Some(0),
        };
    }
    if a.is_empty() || b.is_empty() {
        return rejected("empty_reference");
    }
    // Start narrow: orthologous proteins usually need few indels. The cost
    // bound below expands the band when needed, without changing exactness.
    let mut width = a.len().abs_diff(b.len()) + 16;
    let Some((mut f, mut r, mut cost)) = align(a, b, width) else {
        return rejected("alignment_resource_limit");
    };
    // Any path no worse than this one has at most cost/2 gap residues, so
    // cannot leave this band. Expand until the result is globally exact.
    if width < (cost / 2) as usize {
        width = (cost / 2) as usize;
        let Some(result) = align(a, b, width) else {
            return rejected("alignment_resource_limit");
        };
        (f, r, cost) = result;
    }
    let mut positions = vec![None; a.len()];
    for i in 0..a.len() {
        let mut unique = None;
        let mut alternatives = false;
        for j in f.start[i]..f.end[i] {
            let prefix = f.get(i, j);
            if prefix + 2 + r.get(i + 1, j) == cost {
                alternatives = true;
            }
            if j < b.len() && prefix + u32::from(a[i] != b[j]) + r.get(i + 1, j + 1) == cost {
                if unique.is_some_and(|old| old != j) {
                    alternatives = true;
                }
                unique = Some(j);
            }
        }
        if !alternatives {
            positions[i] = unique;
        }
    }
    let paired = positions.iter().flatten().count();
    let identical = positions
        .iter()
        .enumerate()
        .filter(|(i, j)| j.is_some_and(|j| a[*i] == b[j] && !b"X*".contains(&a[*i])))
        .count();
    if paired == 0 || identical * 5 < paired * 3 {
        return rejected("weak_reference_homology");
    }
    let original = positions.clone();
    for (i, slot) in positions.iter_mut().enumerate() {
        if slot.is_none() {
            continue;
        }
        let start = i.saturating_sub(3);
        let end = (i + 4).min(a.len());
        let matches = (start..end)
            .filter(|&k| original[k].is_some_and(|j| a[k] == b[j] && !b"X*".contains(&a[k])))
            .count();
        if b"X*".contains(&a[i]) || matches * 10 < (end - start) * 7 {
            *slot = None;
        }
    }
    ResidueMap {
        positions,
        status: "unambiguous_anchored_alignment",
        cost: Some(cost),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cost_certified_band_agrees_with_full_matrix_on_large_internal_indels() {
        let source = "ACDEFGHIKLMNPQRSTVWY".repeat(10);
        let target = format!(
            "{}{}{}{}",
            &source[..30],
            &source[70..150],
            "W".repeat(40),
            &source[150..]
        );
        let a = source.as_bytes();
        let b = target.as_bytes();
        let (_, _, upper) = align(a, b, 16).unwrap();
        let width = 16.max((upper / 2) as usize);
        let (f, r, cost) = align(a, b, width).unwrap();
        let (full_f, full_r, full_cost) = align(a, b, a.len().max(b.len())).unwrap();
        assert_eq!(cost, full_cost);
        for (i, base) in a.iter().enumerate() {
            let edges = |forward: &Band, backward: &Band| {
                let mut alternatives = Vec::new();
                for j in 0..=b.len() {
                    if forward.get(i, j) + 2 + backward.get(i + 1, j) == cost {
                        alternatives.push(None);
                    }
                    if j < b.len()
                        && forward.get(i, j) + u32::from(*base != b[j]) + backward.get(i + 1, j + 1)
                            == cost
                    {
                        alternatives.push(Some(j));
                    }
                }
                alternatives
            };
            assert_eq!(edges(&f, &r), edges(&full_f, &full_r));
        }
    }
    #[test]
    fn insertion_preserves_downstream_coordinates_and_is_not_copied() {
        let result = map("MKWVTFGQISLLFLFSSAYSR", "MKWVTFISLLFLFSSAYSR");
        assert_eq!(result.positions[6], None);
        assert_eq!(result.positions[7], None);
        assert_eq!(*result.positions.last().unwrap(), Some(18));
    }
    #[test]
    fn repeat_shift_alternatives_are_unmapped() {
        let result = map("MKWVAAAAAGTFISP", "MKWVAAAAGTFISP");
        assert!(result.positions[4..9].iter().all(Option::is_none));
        assert_eq!(*result.positions.last().unwrap(), Some(13));
    }
    #[test]
    fn substitutions_can_map_when_flanks_agree_but_unrelated_sequences_cannot() {
        let a = "MKWVTFISLLFLFSSAYSR";
        let b = "MKWVTFISLAFLFSSAYSR";
        assert_eq!(map(a, b).positions[9], Some(9));
        assert!(map("AAAAAAAAAAA", "WWWWWWWWWWW")
            .positions
            .iter()
            .all(Option::is_none));
        assert_eq!(
            map(a, a).positions,
            (0..a.len()).map(Some).collect::<Vec<_>>()
        );
    }
}
