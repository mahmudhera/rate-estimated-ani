//! Positional FracMinHash sketches.
//!
//! A [`Sketch`] keeps every sampled k-mer occurrence of a genome in two views:
//!
//! - **Scan view**, in position order: per entry, the hash of the canonical k-mer, its
//!   0-based start within its contig, and whether the forward k-mer is the canonical one.
//!   Each contig's entries are contiguous; `contig_offsets` marks the boundaries.
//! - **Set view**: `by_hash`, a permutation of entry indices sorted by `(hash, index)`.
//!   Equal hashes (repeats) form contiguous runs. Comparisons merge-join two set views.
//!
//! Repeated k-mers are kept: every occurrence is an entry. Masking is a comparison-time choice.

pub mod format;
pub mod kmer;
#[cfg(target_arch = "x86_64")]
pub mod simd;

use std::ops::Range;

use anyhow::{Result, bail, ensure};

use crate::params::{MAX_K, SketchParams};
use kmer::{fmh_threshold, for_each_sampled_kmer};

/// A fixed-length vector of bits, packed into `u64` words (least significant bit first).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BitVec {
    words: Vec<u64>,
    len: usize,
}

impl BitVec {
    /// Rebuilds a bit vector from its words. Fails if the length doesn't match the word count
    /// or bits beyond `len` are set.
    pub fn from_words(words: Vec<u64>, len: usize) -> Result<Self> {
        ensure!(
            words.len() == len.div_ceil(64),
            "bit vector of {len} bits needs {} words, got {}",
            len.div_ceil(64),
            words.len()
        );
        if len % 64 != 0 {
            let last = *words.last().expect("non-empty when len % 64 != 0");
            ensure!(
                last >> (len % 64) == 0,
                "bits set beyond the end of the bit vector"
            );
        }
        Ok(Self { words, len })
    }

    /// Appends a bit.
    pub fn push(&mut self, bit: bool) {
        if self.len % 64 == 0 {
            self.words.push(0);
        }
        if bit {
            *self.words.last_mut().expect("pushed above") |= 1 << (self.len % 64);
        }
        self.len += 1;
    }

    /// The bit at `i`. Panics if `i` is out of range.
    #[inline]
    pub fn get(&self, i: usize) -> bool {
        assert!(
            i < self.len,
            "bit index {i} out of range for length {}",
            self.len
        );
        self.words[i / 64] >> (i % 64) & 1 == 1
    }

    /// Number of bits.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether there are no bits.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The packed words.
    pub fn words(&self) -> &[u64] {
        &self.words
    }
}

/// A positional FracMinHash sketch of one genome. See the module docs for the layout.
#[derive(Clone, Debug, PartialEq)]
pub struct Sketch {
    name: String,
    params: SketchParams,
    contig_names: Vec<String>,
    contig_lengths: Vec<u32>,
    /// `contig_offsets[i]..contig_offsets[i + 1]` are contig `i`'s entries; length is
    /// `num_contigs + 1`.
    contig_offsets: Vec<u32>,
    hashes: Vec<u64>,
    positions: Vec<u32>,
    forward: BitVec,
    by_hash: Vec<u32>,
    skipped_contigs: u32,
    skipped_bases: u64,
}

impl Sketch {
    /// Name of the sketched genome: the input path as given.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Parameters the sketch was built with.
    pub fn params(&self) -> &SketchParams {
        &self.params
    }

    /// Number of contigs kept (at least `min_contig_len` bases long).
    pub fn num_contigs(&self) -> usize {
        self.contig_names.len()
    }

    /// Names (FASTA headers) of the kept contigs.
    pub fn contig_names(&self) -> &[String] {
        &self.contig_names
    }

    /// Lengths of the kept contigs, in bases.
    pub fn contig_lengths(&self) -> &[u32] {
        &self.contig_lengths
    }

    /// Total length of the kept contigs. This is the denominator for aligned fraction.
    pub fn total_length(&self) -> u64 {
        self.contig_lengths.iter().map(|&l| u64::from(l)).sum()
    }

    /// Number of contigs left out for being shorter than `min_contig_len`.
    pub fn skipped_contigs(&self) -> u32 {
        self.skipped_contigs
    }

    /// Total length of the contigs left out.
    pub fn skipped_bases(&self) -> u64 {
        self.skipped_bases
    }

    /// Number of sampled k-mer occurrences.
    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    /// Whether no k-mer was sampled.
    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

    /// Hashes in position order (scan view).
    pub fn hashes(&self) -> &[u64] {
        &self.hashes
    }

    /// 0-based k-mer start positions within their contig, in position order (scan view).
    pub fn positions(&self) -> &[u32] {
        &self.positions
    }

    /// Whether entry `i`'s forward k-mer is its canonical k-mer.
    pub fn is_forward(&self, i: usize) -> bool {
        self.forward.get(i)
    }

    /// Strand bits for all entries.
    pub fn forward_bits(&self) -> &BitVec {
        &self.forward
    }

    /// Entry indices belonging to contig `contig`.
    pub fn contig_entries(&self, contig: usize) -> Range<usize> {
        self.contig_offsets[contig] as usize..self.contig_offsets[contig + 1] as usize
    }

    /// Boundaries of each contig's entries (length `num_contigs + 1`).
    pub fn contig_offsets(&self) -> &[u32] {
        &self.contig_offsets
    }

    /// Entry indices sorted by `(hash, index)` (set view).
    pub fn by_hash(&self) -> &[u32] {
        &self.by_hash
    }

    /// Checks every structural invariant. Sketches read from disk are validated with this.
    pub fn validate(&self) -> Result<()> {
        let k = self.params.k;
        ensure!((1..=MAX_K).contains(&k), "k = {k} is out of range");
        ensure!(self.params.c >= 1, "scale c must be at least 1");
        let n = self.hashes.len();
        ensure!(
            self.positions.len() == n && self.forward.len() == n && self.by_hash.len() == n,
            "entry arrays have different lengths"
        );
        ensure!(n <= u32::MAX as usize, "too many entries");
        let num_contigs = self.contig_names.len();
        ensure!(
            self.contig_lengths.len() == num_contigs
                && self.contig_offsets.len() == num_contigs + 1,
            "contig arrays have inconsistent lengths"
        );
        ensure!(
            self.contig_offsets[0] == 0 && self.contig_offsets[num_contigs] as usize == n,
            "contig offsets don't cover the entries"
        );
        let threshold = fmh_threshold(self.params.c);
        for ci in 0..num_contigs {
            let (start, end) = (self.contig_offsets[ci], self.contig_offsets[ci + 1]);
            ensure!(start <= end, "contig offsets are not monotone");
            let len = u64::from(self.contig_lengths[ci]);
            let mut prev: Option<u32> = None;
            for i in start as usize..end as usize {
                let pos = self.positions[i];
                ensure!(
                    u64::from(pos) + u64::from(k) <= len,
                    "k-mer at {pos} runs past the end of contig {ci}"
                );
                ensure!(
                    prev.is_none_or(|p| p < pos),
                    "positions in contig {ci} are not strictly increasing"
                );
                prev = Some(pos);
                ensure!(
                    self.hashes[i] < threshold,
                    "hash above the FracMinHash threshold"
                );
            }
        }
        // by_hash must be a permutation sorted by (hash, index).
        let mut seen = vec![false; n];
        for &i in &self.by_hash {
            let i = i as usize;
            ensure!(
                i < n && !std::mem::replace(&mut seen[i], true),
                "set view is not a permutation"
            );
        }
        ensure!(
            self.by_hash
                .windows(2)
                .all(|w| (self.hashes[w[0] as usize], w[0]) < (self.hashes[w[1] as usize], w[1])),
            "set view is not sorted by (hash, index)"
        );
        Ok(())
    }
}

/// Builds a [`Sketch`] one contig at a time.
#[derive(Debug)]
pub struct SketchBuilder {
    name: String,
    params: SketchParams,
    threshold: u64,
    contig_names: Vec<String>,
    contig_lengths: Vec<u32>,
    contig_offsets: Vec<u32>,
    hashes: Vec<u64>,
    positions: Vec<u32>,
    forward: BitVec,
    skipped_contigs: u32,
    skipped_bases: u64,
}

impl SketchBuilder {
    /// Starts an empty sketch named `name`.
    pub fn new(name: impl Into<String>, params: SketchParams) -> Self {
        assert!((1..=MAX_K).contains(&params.k), "k must be in 1..={MAX_K}");
        Self {
            name: name.into(),
            params,
            threshold: fmh_threshold(params.c),
            contig_names: Vec::new(),
            contig_lengths: Vec::new(),
            contig_offsets: vec![0],
            hashes: Vec::new(),
            positions: Vec::new(),
            forward: BitVec::default(),
            skipped_contigs: 0,
            skipped_bases: 0,
        }
    }

    /// Sketches one contig. Contigs shorter than `min_contig_len` are counted as skipped.
    pub fn add_contig(&mut self, name: &str, seq: &[u8]) -> Result<()> {
        let len = seq.len();
        if len < self.params.min_contig_len as usize {
            self.skipped_contigs += 1;
            self.skipped_bases += len as u64;
            return Ok(());
        }
        let Ok(len32) = u32::try_from(len) else {
            bail!(
                "contig {name} is {len} bases long; the limit is {}",
                u32::MAX
            );
        };
        self.contig_names.push(name.to_owned());
        self.contig_lengths.push(len32);
        // Positions fit in u32 because the contig does.
        for_each_sampled_kmer(seq, self.params.k, self.threshold, |s| {
            self.hashes.push(s.hash);
            self.positions.push(s.pos as u32);
            self.forward.push(s.forward);
        });
        let Ok(end) = u32::try_from(self.hashes.len()) else {
            bail!("more than {} sampled k-mers in one genome", u32::MAX);
        };
        self.contig_offsets.push(end);
        Ok(())
    }

    /// Finishes the sketch and builds the set view.
    pub fn finish(self) -> Sketch {
        // Sort packed (hash, index) keys: contiguous and compared in one instruction, which is
        // about twice as fast as sorting indices by an indirect key.
        let mut keyed: Vec<u128> = self
            .hashes
            .iter()
            .enumerate()
            .map(|(i, &h)| (u128::from(h) << 32) | i as u128)
            .collect();
        keyed.sort_unstable();
        let by_hash: Vec<u32> = keyed.into_iter().map(|x| x as u32).collect();
        Sketch {
            name: self.name,
            params: self.params,
            contig_names: self.contig_names,
            contig_lengths: self.contig_lengths,
            contig_offsets: self.contig_offsets,
            hashes: self.hashes,
            positions: self.positions,
            forward: self.forward,
            by_hash,
            skipped_contigs: self.skipped_contigs,
            skipped_bases: self.skipped_bases,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::kmer::tests::{random_seq, revcomp};
    use super::*;

    fn params(k: u8, c: u32, min_contig_len: u32) -> SketchParams {
        SketchParams {
            k,
            c,
            min_contig_len,
        }
    }

    fn build(contigs: &[(&str, &[u8])], p: SketchParams) -> Sketch {
        let mut b = SketchBuilder::new("test", p);
        for (name, seq) in contigs {
            b.add_contig(name, seq).unwrap();
        }
        b.finish()
    }

    #[test]
    fn bitvec_round_trips() {
        let mut b = BitVec::default();
        let bits: Vec<bool> = (0..130).map(|i| i % 3 == 0).collect();
        for &x in &bits {
            b.push(x);
        }
        assert_eq!(b.len(), 130);
        for (i, &x) in bits.iter().enumerate() {
            assert_eq!(b.get(i), x);
        }
        let again = BitVec::from_words(b.words().to_vec(), b.len()).unwrap();
        assert_eq!(again, b);
        assert!(BitVec::from_words(vec![u64::MAX], 3).is_err());
        assert!(BitVec::from_words(vec![0, 0], 3).is_err());
    }

    #[test]
    fn layout_invariants_hold() {
        let a = random_seq(20_000, 3);
        let b = random_seq(9_000, 4);
        let s = build(
            &[("a", &a), ("tiny", b"ACGT"), ("b", &b)],
            params(15, 20, 100),
        );
        s.validate().unwrap();
        assert_eq!(s.num_contigs(), 2);
        assert_eq!(s.contig_names(), ["a", "b"]);
        assert_eq!(s.total_length(), 29_000);
        assert_eq!((s.skipped_contigs(), s.skipped_bases()), (1, 4));
        assert_eq!(s.contig_entries(0).end, s.contig_entries(1).start);
        assert_eq!(s.contig_entries(1).end, s.len());
        // Positions restart at each contig.
        assert!(s.positions()[s.contig_entries(1).start] < 9_000);
    }

    #[test]
    fn repeats_are_kept_and_grouped_in_the_set_view() {
        // Three copies separated by N, so no k-mer spans two copies.
        let unit = random_seq(1_000, 5);
        let mut genome = unit.clone();
        for _ in 0..2 {
            genome.push(b'N');
            genome.extend_from_slice(&unit);
        }
        let s = build(&[("rep", &genome)], params(17, 10, 0));
        s.validate().unwrap();
        // Every hash occurs exactly three times, and the set view groups the copies.
        for run in s.by_hash().chunks(3) {
            let h = s.hashes()[run[0] as usize];
            assert!(run.iter().all(|&i| s.hashes()[i as usize] == h));
        }
        assert_eq!(s.len() % 3, 0);
    }

    #[test]
    fn reverse_complement_gives_mirrored_sketch() {
        let seq = random_seq(50_000, 11);
        let rc = revcomp(&seq);
        let p = params(17, 25, 0);
        let f = build(&[("f", &seq)], p);
        let r = build(&[("r", &rc)], p);
        assert_eq!(f.len(), r.len());
        let mut hf = f.hashes().to_vec();
        let mut hr = r.hashes().to_vec();
        hf.sort_unstable();
        hr.sort_unstable();
        assert_eq!(hf, hr);
        // Entry i on the forward strand is entry n-1-i on the reverse complement.
        let n = f.len();
        for i in 0..n {
            let j = n - 1 - i;
            assert_eq!(f.hashes()[i], r.hashes()[j]);
            assert_eq!(
                f.positions()[i] as usize,
                seq.len() - 17 - r.positions()[j] as usize
            );
            assert_eq!(f.is_forward(i), !r.is_forward(j));
        }
    }

    #[test]
    fn empty_genome_is_valid() {
        let s = build(&[], SketchParams::default());
        s.validate().unwrap();
        assert!(s.is_empty());
        assert_eq!(s.total_length(), 0);
    }

    #[test]
    fn validate_catches_broken_invariants() {
        let seq = random_seq(10_000, 8);
        let good = build(&[("a", &seq)], params(15, 10, 0));
        let mut bad = good.clone();
        bad.by_hash.swap(0, 1);
        assert!(bad.validate().is_err());
        let mut bad = good.clone();
        bad.positions[1] = bad.positions[0];
        assert!(bad.validate().is_err());
        let mut bad = good.clone();
        bad.contig_offsets[1] -= 1;
        assert!(bad.validate().is_err());
        let mut bad = good;
        bad.hashes[0] = u64::MAX;
        assert!(bad.validate().is_err());
    }
}
