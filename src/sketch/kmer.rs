//! k-mer extraction and FracMinHash sampling.
//!
//! Adapted from skani v0.4.0 (`src/seeding.rs`, `src/types.rs`; MIT OR Apache-2.0,
//! Copyright (c) 2022 Jim Shaw). The hash is `mm_hash64`, the invertible 64-bit integer hash
//! that skani takes from minimap2 (MIT, Heng Li). Unlike skani, the k-mer hashed here is
//! always the canonical k-mer at the reported position, and any non-ACGT base ends a k-mer.

use crate::params::MAX_K;

/// 2-bit code of each byte: `A/a` 0, `C/c` 1, `G/g` 2, `T/t` 3. Every other byte maps to 4,
/// which means "not a base" and breaks k-mers.
pub const NT4: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// Invertible 64-bit integer hash (Thomas Wang's, as used by minimap2 and skani).
#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

/// FracMinHash keeps a k-mer when its hash is strictly below this value, which keeps about
/// one in `c` k-mers. Same rule as skani.
#[inline]
pub fn fmh_threshold(c: u32) -> u64 {
    assert!(c >= 1, "FracMinHash scale must be at least 1");
    u64::MAX / u64::from(c)
}

/// One sampled k-mer occurrence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SampledKmer {
    /// 0-based start of the k-mer in its sequence.
    pub pos: usize,
    /// Hash of the canonical k-mer.
    pub hash: u64,
    /// True when the k-mer as read on the forward strand is the canonical one.
    pub forward: bool,
}

/// Calls `emit` for every k-mer of `seq` whose canonical hash is below `threshold`, in
/// increasing position order. k-mers that contain a non-ACGT byte are skipped.
///
/// The canonical k-mer is the smaller of the forward k-mer and its reverse complement, in
/// 2-bit encoding with the first base in the most significant bits.
///
/// Uses the AVX2 path when the CPU supports it, unless the environment variable
/// `KADENCE_NO_SIMD` is set; both paths give identical output.
pub fn for_each_sampled_kmer<F: FnMut(SampledKmer)>(
    seq: &[u8],
    k: u8,
    threshold: u64,
    mut emit: F,
) {
    #[cfg(target_arch = "x86_64")]
    if use_simd() {
        assert!((1..=MAX_K).contains(&k), "k must be in 1..={MAX_K}");
        let k = usize::from(k);
        let mut held: [Vec<SampledKmer>; 3] = Default::default();
        let mut i = 0;
        // Split into maximal runs of ACGT; k-mers never span a run boundary.
        while i < seq.len() {
            while i < seq.len() && NT4[usize::from(seq[i])] > 3 {
                i += 1;
            }
            let start = i;
            while i < seq.len() && NT4[usize::from(seq[i])] <= 3 {
                i += 1;
            }
            if i - start >= k {
                // SAFETY: use_simd() checked that the CPU supports AVX2.
                #[allow(unsafe_code)]
                unsafe {
                    super::simd::sample_run(
                        &seq[start..i],
                        start,
                        k,
                        threshold,
                        &mut held,
                        &mut emit,
                    );
                }
            }
        }
        return;
    }
    for_each_sampled_kmer_scalar(seq, k, threshold, emit);
}

/// Whether to use the AVX2 path: supported by the CPU and not disabled by `KADENCE_NO_SIMD`.
#[cfg(target_arch = "x86_64")]
fn use_simd() -> bool {
    static USE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *USE.get_or_init(|| super::simd::available() && std::env::var_os("KADENCE_NO_SIMD").is_none())
}

/// Portable scalar implementation of [`for_each_sampled_kmer`]. It is the fallback on CPUs
/// without AVX2 and the reference the AVX2 path is tested against.
pub fn for_each_sampled_kmer_scalar<F: FnMut(SampledKmer)>(
    seq: &[u8],
    k: u8,
    threshold: u64,
    mut emit: F,
) {
    assert!((1..=MAX_K).contains(&k), "k must be in 1..={MAX_K}");
    let k = usize::from(k);
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let rc_shift = 2 * (k - 1);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    // Number of consecutive valid bases ending at the current position.
    let mut run: usize = 0;
    for (i, &byte) in seq.iter().enumerate() {
        let code = NT4[usize::from(byte)];
        if code > 3 {
            run = 0;
            continue;
        }
        let code = u64::from(code);
        fwd = ((fwd << 2) | code) & mask;
        rev = (rev >> 2) | ((3 ^ code) << rc_shift);
        run += 1;
        if run >= k {
            let (canonical, forward) = if fwd <= rev {
                (fwd, true)
            } else {
                (rev, false)
            };
            let hash = mm_hash64(canonical);
            if hash < threshold {
                emit(SampledKmer {
                    pos: i + 1 - k,
                    hash,
                    forward,
                });
            }
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Reverse complement of an ACGT string (other bytes become `N`).
    pub(crate) fn revcomp(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' | b'a' => b'T',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'T' | b't' => b'A',
                _ => b'N',
            })
            .collect()
    }

    fn encode(kmer: &[u8]) -> Option<u64> {
        kmer.iter()
            .try_fold(0u64, |acc, &b| match NT4[usize::from(b)] {
                c @ 0..=3 => Some((acc << 2) | u64::from(c)),
                _ => None,
            })
    }

    /// Brute-force reference: every window, canonical via an explicit reverse complement.
    pub(crate) fn naive_sampled(seq: &[u8], k: u8, threshold: u64) -> Vec<SampledKmer> {
        let k = usize::from(k);
        let mut out = Vec::new();
        if seq.len() < k {
            return out;
        }
        for pos in 0..=seq.len() - k {
            let window = &seq[pos..pos + k];
            let (Some(f), Some(r)) = (encode(window), encode(&revcomp(window))) else {
                continue;
            };
            let (canonical, forward) = if f <= r { (f, true) } else { (r, false) };
            let hash = mm_hash64(canonical);
            if hash < threshold {
                out.push(SampledKmer { pos, hash, forward });
            }
        }
        out
    }

    pub(crate) fn sampled(seq: &[u8], k: u8, threshold: u64) -> Vec<SampledKmer> {
        let mut out = Vec::new();
        for_each_sampled_kmer(seq, k, threshold, |s| out.push(s));
        let mut scalar = Vec::new();
        for_each_sampled_kmer_scalar(seq, k, threshold, |s| scalar.push(s));
        assert_eq!(out, scalar, "dispatched and scalar paths differ");
        out
    }

    /// Deterministic pseudo-random ACGT sequence (xorshift), so tests need no extra crates.
    pub(crate) fn random_seq(len: usize, mut state: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize & 3]
            })
            .collect()
    }

    #[test]
    fn encoding_table() {
        for (b, c) in [(b'A', 0), (b'c', 1), (b'G', 2), (b't', 3)] {
            assert_eq!(NT4[usize::from(b)], c);
        }
        for b in *b"NnRYU-\n " {
            assert_eq!(NT4[usize::from(b)], 4, "{}", b as char);
        }
    }

    #[test]
    fn hash_is_a_bijection_on_samples() {
        // Distinct inputs never collide (mm_hash64 is invertible).
        let mut seen = std::collections::HashSet::new();
        for x in 0..100_000u64 {
            assert!(seen.insert(mm_hash64(x)));
        }
    }

    #[test]
    fn rolling_matches_brute_force() {
        let mut seq = random_seq(5_000, 42);
        // Sprinkle in non-ACGT bytes and lowercase bases.
        for (i, b) in [
            (10, b'N'),
            (11, b'N'),
            (700, b'R'),
            (2_500, b'n'),
            (4_999, b'Y'),
        ] {
            seq[i] = b;
        }
        for b in &mut seq[1_000..1_100] {
            *b = b.to_ascii_lowercase();
        }
        for k in [1, 2, 5, 15, 16, 17, 21, 31] {
            for c in [1, 7, 125] {
                let t = fmh_threshold(c);
                assert_eq!(
                    sampled(&seq, k, t),
                    naive_sampled(&seq, k, t),
                    "k={k} c={c}"
                );
            }
        }
    }

    #[test]
    fn non_acgt_bases_break_kmers() {
        let k = 5;
        let seq = b"ACGTACGTRACGTACGT";
        let all = sampled(seq, k, u64::MAX);
        assert!(all.iter().all(|s| !(s.pos..s.pos + 5).contains(&8)));
        // 4 k-mers fit on each side of the R.
        assert_eq!(all.len(), 8);
    }

    #[test]
    fn short_sequences_yield_nothing() {
        assert!(sampled(b"ACG", 5, u64::MAX).is_empty());
        assert!(sampled(b"", 5, u64::MAX).is_empty());
        assert_eq!(sampled(b"ACGTA", 5, u64::MAX).len(), 1);
    }

    #[test]
    fn canonical_hashes_are_strand_symmetric() {
        let seq = random_seq(3_000, 9);
        let rc = revcomp(&seq);
        let k = 17;
        let f = sampled(&seq, k, fmh_threshold(5));
        let r = sampled(&rc, k, fmh_threshold(5));
        assert_eq!(f.len(), r.len());
        // The occurrence at forward position p appears at p' = L - k - p on the reverse
        // complement, with the same hash and the opposite orientation (odd k: no palindromes).
        let len = seq.len();
        let mut mirrored: Vec<SampledKmer> = r
            .iter()
            .map(|s| SampledKmer {
                pos: len - usize::from(k) - s.pos,
                hash: s.hash,
                forward: !s.forward,
            })
            .collect();
        mirrored.sort_by_key(|s| s.pos);
        assert_eq!(f, mirrored);
    }

    #[test]
    fn sampling_rate_is_one_in_c() {
        let seq = random_seq(2_000_000, 1);
        let c = 125;
        let n = sampled(&seq, 17, fmh_threshold(c)).len() as f64;
        let expected = (seq.len() - 16) as f64 / f64::from(c);
        // Binomial sd is about sqrt(16_000) = 126; allow 5 sd.
        assert!(
            (n - expected).abs() < 5.0 * expected.sqrt(),
            "{n} vs {expected}"
        );
    }
}
