//! AVX2 k-mer sampling, selected at run time on x86_64 CPUs that support it.
//!
//! Ported from skani v0.4.0 (`src/avx2_seeding.rs`; MIT OR Apache-2.0, Copyright (c) 2022
//! Jim Shaw) to 64-bit hashes and k up to 31. As in skani, an ACGT run is split into four
//! contiguous lanes that roll and hash in parallel. Unlike skani, each lane's canonical
//! k-mer is the true canonical k-mer at its position, exactly as in the scalar path
//! ([`super::kmer::for_each_sampled_kmer`]); a test checks that the two paths agree.

// Intrinsics need `unsafe` on older compilers but not on newer ones.
#![allow(unsafe_code, unused_unsafe)]

use std::arch::x86_64::*;

use super::kmer::{NT4, SampledKmer, mm_hash64};

/// Whether the AVX2 path can run on this CPU.
pub fn available() -> bool {
    is_x86_feature_detected!("avx2")
}

/// AVX2 version of `mm_hash64` on four lanes.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn mm_hash256(kmer: __m256i) -> __m256i {
    unsafe {
        let mut key = kmer;
        // key = !(key + (key << 21))
        key = _mm256_add_epi64(key, _mm256_slli_epi64::<21>(key));
        key = _mm256_xor_si256(key, _mm256_set1_epi64x(-1));
        key = _mm256_xor_si256(key, _mm256_srli_epi64::<24>(key));
        key = _mm256_add_epi64(
            _mm256_add_epi64(key, _mm256_slli_epi64::<3>(key)),
            _mm256_slli_epi64::<8>(key),
        );
        key = _mm256_xor_si256(key, _mm256_srli_epi64::<14>(key));
        key = _mm256_add_epi64(
            _mm256_add_epi64(key, _mm256_slli_epi64::<2>(key)),
            _mm256_slli_epi64::<4>(key),
        );
        key = _mm256_xor_si256(key, _mm256_srli_epi64::<28>(key));
        _mm256_add_epi64(key, _mm256_slli_epi64::<31>(key))
    }
}

/// Samples the k-mers of `run` (ACGT only; starts at `offset` in its contig), emitting them
/// in position order. `held` is scratch space and must be empty.
///
/// # Safety
///
/// The CPU must support AVX2 (check [`available`]).
#[target_feature(enable = "avx2")]
pub unsafe fn sample_run<F: FnMut(SampledKmer)>(
    run: &[u8],
    offset: usize,
    k: usize,
    threshold: u64,
    held: &mut [Vec<SampledKmer>; 3],
    emit: &mut F,
) {
    debug_assert!(run.len() >= k && held.iter().all(Vec::is_empty));
    let m = run.len() - k + 1; // number of k-mers
    let q = m / 4; // k-mers per lane; the remaining m % 4 go through the scalar tail
    let mask = (1u64 << (2 * k)) - 1;
    let rc_shift = 2 * (k - 1);

    if q > 0 {
        unsafe {
            // Lane l covers k-mers l*q .. (l+1)*q. Prime each lane with its first k - 1 bases.
            let mut fwd = [0u64; 4];
            let mut rev = [0u64; 4];
            for l in 0..4 {
                for &b in &run[l * q..l * q + k - 1] {
                    let code = u64::from(NT4[usize::from(b)]);
                    fwd[l] = (fwd[l] << 2) | code;
                    rev[l] = (rev[l] >> 2) | ((3 ^ code) << rc_shift);
                }
            }
            let mut vf =
                _mm256_set_epi64x(fwd[3] as i64, fwd[2] as i64, fwd[1] as i64, fwd[0] as i64);
            let mut vr =
                _mm256_set_epi64x(rev[3] as i64, rev[2] as i64, rev[1] as i64, rev[0] as i64);
            let vmask = _mm256_set1_epi64x(mask as i64);
            let three = _mm256_set1_epi64x(3);
            let shift = _mm256_set1_epi64x(rc_shift as i64);
            // Unsigned hash < threshold, via signed compare with the sign bits flipped.
            let sign = _mm256_set1_epi64x(i64::MIN);
            let vthr = _mm256_xor_si256(_mm256_set1_epi64x(threshold as i64), sign);
            let mut hashes = [0u64; 4];

            for j in 0..q {
                let at = k - 1 + j;
                let codes = _mm256_set_epi64x(
                    i64::from(NT4[usize::from(run[3 * q + at])]),
                    i64::from(NT4[usize::from(run[2 * q + at])]),
                    i64::from(NT4[usize::from(run[q + at])]),
                    i64::from(NT4[usize::from(run[at])]),
                );
                vf = _mm256_and_si256(_mm256_or_si256(_mm256_slli_epi64::<2>(vf), codes), vmask);
                vr = _mm256_or_si256(
                    _mm256_srli_epi64::<2>(vr),
                    _mm256_sllv_epi64(_mm256_xor_si256(codes, three), shift),
                );
                // k-mers are below 2^62, so the signed compare is exact.
                let fwd_gt = _mm256_cmpgt_epi64(vf, vr);
                let canonical = _mm256_blendv_epi8(vf, vr, fwd_gt);
                let h = mm_hash256(canonical);
                let hit = _mm256_cmpgt_epi64(vthr, _mm256_xor_si256(h, sign));
                let hits = _mm256_movemask_pd(_mm256_castsi256_pd(hit));
                if hits != 0 {
                    _mm256_storeu_si256(hashes.as_mut_ptr().cast(), h);
                    let gt = _mm256_movemask_pd(_mm256_castsi256_pd(fwd_gt));
                    for l in 0..4 {
                        if hits >> l & 1 == 1 {
                            let s = SampledKmer {
                                pos: offset + l * q + j,
                                hash: hashes[l],
                                forward: gt >> l & 1 == 0,
                            };
                            if l == 0 {
                                emit(s);
                            } else {
                                held[l - 1].push(s);
                            }
                        }
                    }
                }
            }
        }
        for lane in held.iter_mut() {
            for s in lane.drain(..) {
                emit(s);
            }
        }
    }

    // Scalar tail: the last m % 4 k-mers (all of them when the run is short).
    let first = 4 * q;
    let (mut fwd, mut rev) = (0u64, 0u64);
    for &b in &run[first..first + k - 1] {
        let code = u64::from(NT4[usize::from(b)]);
        fwd = (fwd << 2) | code;
        rev = (rev >> 2) | ((3 ^ code) << rc_shift);
    }
    for (j, &b) in run[first + k - 1..].iter().enumerate() {
        let code = u64::from(NT4[usize::from(b)]);
        fwd = ((fwd << 2) | code) & mask;
        rev = (rev >> 2) | ((3 ^ code) << rc_shift);
        let hash = mm_hash64(fwd.min(rev));
        if hash < threshold {
            emit(SampledKmer {
                pos: offset + first + j,
                hash,
                forward: fwd <= rev,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::kmer::tests::{naive_sampled, random_seq};
    use super::super::kmer::{fmh_threshold, for_each_sampled_kmer_scalar};
    use super::*;

    fn simd_sampled(seq: &[u8], k: u8, threshold: u64) -> Vec<SampledKmer> {
        let mut out = Vec::new();
        let mut held: [Vec<SampledKmer>; 3] = Default::default();
        let k = usize::from(k);
        // Same run splitting as the dispatcher.
        let mut i = 0;
        while i < seq.len() {
            while i < seq.len() && NT4[usize::from(seq[i])] > 3 {
                i += 1;
            }
            let start = i;
            while i < seq.len() && NT4[usize::from(seq[i])] <= 3 {
                i += 1;
            }
            if i - start >= k {
                unsafe {
                    sample_run(&seq[start..i], start, k, threshold, &mut held, &mut |s| {
                        out.push(s)
                    })
                };
            }
        }
        out
    }

    #[test]
    fn hash256_matches_scalar() {
        if !available() {
            return;
        }
        let xs = [0u64, 1, 12345, u64::MAX >> 2];
        let v = unsafe {
            mm_hash256(_mm256_set_epi64x(
                xs[3] as i64,
                xs[2] as i64,
                xs[1] as i64,
                xs[0] as i64,
            ))
        };
        let mut out = [0u64; 4];
        unsafe { _mm256_storeu_si256(out.as_mut_ptr().cast(), v) };
        assert_eq!(out, xs.map(mm_hash64));
    }

    #[test]
    fn simd_matches_scalar_and_brute_force() {
        if !available() {
            return;
        }
        let mut seq = random_seq(20_000, 77);
        for (i, b) in [
            (5, b'N'),
            (300, b'R'),
            (301, b'N'),
            (9_999, b'n'),
            (19_990, b'Y'),
        ] {
            seq[i] = b;
        }
        for k in [1u8, 2, 3, 15, 16, 17, 21, 31] {
            for c in [1u32, 3, 125] {
                let t = fmh_threshold(c);
                let mut scalar = Vec::new();
                for_each_sampled_kmer_scalar(&seq, k, t, |s| scalar.push(s));
                assert_eq!(simd_sampled(&seq, k, t), scalar, "k={k} c={c}");
            }
        }
        assert_eq!(
            simd_sampled(&seq[..3_000], 17, fmh_threshold(2)),
            naive_sampled(&seq[..3_000], 17, fmh_threshold(2))
        );
    }

    #[test]
    fn short_runs_use_the_tail() {
        if !available() {
            return;
        }
        // Runs with 1..12 k-mers exercise q = 0..3 and every remainder.
        for len in 17..29 {
            let seq = random_seq(len, len as u64);
            assert_eq!(
                simd_sampled(&seq, 17, u64::MAX),
                naive_sampled(&seq, 17, u64::MAX),
                "len={len}"
            );
        }
    }
}
