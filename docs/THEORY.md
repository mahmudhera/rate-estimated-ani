# Theory notes: rate-based ANI and AF from positional FracMinHash sketches

Working name: **kadence** (k-mer cadence). Every number below is reproduced by
[`analysis/theory.py`](../analysis/theory.py) (`python3 analysis/theory.py`), assuming
c = 125, a 5 Mbp genome, and random sequence where stated.

Notation: `a` = ANI, `k` = k-mer size, `c` = FracMinHash scale (keep 1/c of k-mers),
`x` = **patience** (consecutive non-matching sketched k-mers before a region closes),
`y` = **min-anchors** (matches required to accept a region), `q = a^k`, `s` = spurious match rate.

## Verdict

**Promising, no fatal flaw.** The method belongs to the same estimator family as skani:
ANI = (matched / sketched k-mers inside homologous regions)^(1/k). Only the region finder
differs. skani runs banded 2-D chaining, using coordinates on both genomes, inside 20 kb query
chunks. We run 1-D run segmentation along one genome's sketch.

Taken literally, the design has four problems. Each has an O(n) fix:

| # | Issue | Severity | Fix |
|---|---|---|---|
| 1 | Measuring the rate per base pair throws away information | Moderate: 7× variance at 99% ANI, 67× at 99.9% | Measure the rate per sketched k-mer. We store positions, so we know the count. |
| 2 | No single fixed patience works from 80% to 99.9% ANI | Major, for AF | Scale patience to the pair's own match rate (adaptive) |
| 3 | Regions are defined by the same matches we count, which biases ANI upward | Major below about 90% ANI | Truncated-geometric MLE of the rate |
| 4 | Without a colinearity check, spurious and paralogous matches count | Major at k = 15 below about 95% ANI | Use k ≥ 17, or carry B's position and require diagonal consistency |

## 1. Model

FracMinHash keeps a k-mer iff `h(kmer) < 2^64 / c`, which depends only on the k-mer itself.
So if a sampled k-mer of A is unchanged in B, B sampled it too. Under i.i.d. substitutions:

- Inside a homologous region, each entry of A's position-ordered sketch is a match with
  probability `q = a^k`. The entries are approximately independent, because sampled k-mers
  rarely overlap when c ≫ k.
- The number of misses between consecutive matches, G, is Geometric(q): `P(G = g) = (1−q)^g q`.
- Measured in base pairs, matches form a thinned Poisson process with rate `q / c`.
  This is the "Poisson rate".
- Outside homology, matches occur at a spurious rate `s ≈ 1 − exp(−2|B| / 4^k)` for random
  sequence. Real genomes have a higher rate because of skewed composition.

## 2. Use the sketch index as the clock, not base pairs

The per-bp estimator `q̂ = m·c / L` has variance `q / n̄`. The count estimator `q̂ = m / n` has
variance `q(1−q) / n`. Because we store every sketched position, we know n exactly.

| ANI | q (k = 15) | Var(per-bp) / Var(count) |
|---|---|---|
| 0.999 | 0.985 | 67× |
| 0.99 | 0.860 | 7.1× |
| 0.95 | 0.463 | 1.9× |
| 0.90 | 0.206 | 1.3× |

The per-bp rate also assumes uniform sampling density, which repeats, runs of N and
low-complexity sequence violate. The rate idea still holds, but its clock is the sketch index.

## 3. Match density sets the resolution

`q = a^k`, with the mean distance between matches (c/q) in brackets, at c = 125:

| ANI | k = 15 | k = 16 | k = 17 | k = 19 | k = 21 |
|---|---|---|---|---|---|
| 0.99 | 0.860 [0.15 kb] | 0.852 [0.15 kb] | 0.843 [0.15 kb] | 0.826 [0.15 kb] | 0.810 [0.15 kb] |
| 0.95 | 0.463 [0.27 kb] | 0.440 [0.28 kb] | 0.418 [0.30 kb] | 0.377 [0.33 kb] | 0.341 [0.37 kb] |
| 0.90 | 0.206 [0.61 kb] | 0.185 [0.67 kb] | 0.167 [0.75 kb] | 0.135 [0.93 kb] | 0.109 [1.14 kb] |
| 0.85 | 0.087 [1.43 kb] | 0.074 [1.68 kb] | 0.063 [1.98 kb] | 0.046 [2.74 kb] | 0.033 [3.79 kb] |
| 0.80 | 0.035 [3.55 kb] | 0.028 [4.44 kb] | 0.023 [5.55 kb] | 0.014 [8.67 kb] | 0.009 [13.6 kb] |

At c = 125, region boundaries can be resolved at gene scale (about 1 kb) only above about
90% ANI. Every k-mer method at this sampling density has a floor around 80–85% ANI; skani
reports reliable ANI at about 82% and above. Resolution improves linearly as c decreases.

## 4. Patience (x)

A false break happens when a run of misses inside true homology reaches x. The fraction of a
homologous region lost to false breaks is `(1−q)^x · (q·x + 1 − q)`. The table gives the
smallest x that keeps this loss at or below 1%. In base pairs, x·c is also the smallest
non-homologous island the method can detect.

| ANI | k = 15 | k = 16 | k = 17 | k = 19 | k = 21 |
|---|---|---|---|---|---|
| 0.99 | 3 (0.4 kb) | 3 (0.4 kb) | 4 (0.5 kb) | 4 (0.5 kb) | 4 (0.5 kb) |
| 0.95 | 11 (1.4 kb) | 11 (1.4 kb) | 12 (1.5 kb) | 14 (1.8 kb) | 16 (2.0 kb) |
| 0.90 | 29 (3.6 kb) | 32 (4.0 kb) | 36 (4.5 kb) | 46 (5.8 kb) | 57 (7.1 kb) |
| 0.85 | 72 (9.0 kb) | 86 (10.8 kb) | 102 (12.8 kb) | 142 (17.8 kb) | 198 (24.8 kb) |
| 0.80 | 185 (23 kb) | 232 (29 kb) | 291 (36 kb) | 457 (57 kb) | 716 (90 kb) |

**This creates a dilemma.** At k = 17, an 85% ANI pair needs x = 100. Applied to a 99% ANI
pair, x = 100 bridges every non-homologous island shorter than about 12.8 kb and counts it as
aligned. AF is then inflated in strain-level comparisons, where AF matters most. Going the other
way, x = 4 is right at 99% ANI but breaks an 85% pair into fragments.

**Fix: adaptive patience** `x(q) = ⌈ln α / ln(1 − q)⌉`, where α is the per-gap false-break
probability. With α = 10⁻³:

| k | 0.99 | 0.95 | 0.90 | 0.85 | 0.80 |
|---|---|---|---|---|---|
| 15 | 4 | 12 | 30 | 76 | 193 |
| 17 | 4 | 13 | 38 | 106 | 304 |
| 19 | 4 | 15 | 48 | 149 | 476 |

Adaptive patience runs in stages:

1. Segment with a floor patience `x(q_floor)` taken from `--min-ani`.
2. Estimate q̂ as in §5.
3. Re-segment with `x(q̂)`.
4. Repeat steps 2–3 until x stops changing, usually within 2–3 rounds.

Each round is one linear scan over a bit vector, so the cost is negligible.

**Interpretation.** This rule is exactly the Viterbi decision of a two-state HMM with states
"homologous" (match probability q) and "background" (s ≈ 0), and state-switch probability α.
Leaving the homologous state costs `ln(1/α)`, and each miss inside it costs `ln(1/(1−q))`.
Patience is therefore a principled quantity, not just a knob. A full HMM is a natural M5
upgrade: it would handle local variation in ANI and give posterior region boundaries.

## 5. Estimator bias from the segmentation rule

Each region starts and ends with a match, and by construction every internal gap is shorter
than x. So the naive estimate (m/n inside regions)^(1/k) is biased upward unless x ≫ 1/q.

Theory, pooling gaps across regions at k = 17:

| x | 0.95 | 0.90 | 0.85 | 0.80 |
|---|---|---|---|---|
| 10 | 0.9511 | 0.9207 | 0.9100 | 0.9064 |
| 20 | 0.9500 | 0.9050 | 0.8824 | 0.8747 |
| 40 | 0.9500 | 0.9002 | 0.8613 | 0.8450 |
| 80 | 0.9500 | 0.9000 | 0.8514 | 0.8203 |

Monte Carlo with 400k sketched k-mers at k = 17:

| True ANI | x | Naive m/n | Truncated-geometric MLE |
|---|---|---|---|
| 0.85 | 20 | 0.8969 | **0.8505** |
| 0.85 | 40 | 0.8649 | **0.8499** |
| 0.90 | 15 | 0.9128 | **0.8997** |
| 0.95 | 6 | 0.9571 | **0.9499** |

The naive estimate in the Monte Carlo is worse than the pooled-gap theory predicts. Short
regions also contribute a "free" boundary match each, which the pooled-gap formula ignores.

**Fix.** Treat the m−1 internal gaps of each accepted region as draws from Geometric(q),
truncated to G < x. The log-likelihood of one gap is

`ℓ(q; g) = g·ln(1−q) + ln q − ln(1 − (1−q)^x)`.

Setting the score to zero gives `E_q[G | G < x] = Ḡ`, the mean observed gap, where

`E_q[G | G < x] = (1−q)·(1 − x(1−q)^(x−1) + (x−1)(1−q)^x) / (q·(1 − (1−q)^x))`.

This is monotone in q, so bisection solves it. The Fisher information gives an analytic
standard error, and the delta method turns it into an ANI confidence interval. No bootstrap is
needed; skani uses one. This MLE is the cleanest form of the "Poisson rate" idea: the rate is
estimated from the inter-arrival times of matches.

## 6. Spurious matches without colinearity

The spurious rate s per sketched k-mer (random 5 Mbp genome, both strands):

| k | s |
|---|---|
| 15 | 9.3e-3 |
| 16 | 2.3e-3 |
| 17 | 5.8e-4 |
| 19 | 3.6e-5 |
| 21 | 2.3e-6 |

Spurious-region pressure is s·x, with x tuned to the ANI floor. It should be well below 0.1.

| ANI floor | k = 15 | k = 16 | k = 17 | k = 19 | k = 21 |
|---|---|---|---|---|---|
| 0.95 | 0.102 | 0.026 | 0.007 | 0.001 | 0.000 |
| 0.90 | 0.269 | 0.074 | 0.021 | 0.002 | 0.000 |
| 0.85 | 0.667 | 0.200 | 0.059 | 0.005 | 0.000 |
| 0.80 | 1.715 | 0.540 | 0.169 | 0.017 | 0.002 |

Expected false regions per 5 Mbp genome, `≈ (G/c)·s·(s·x)^(y−1)`, with x tuned for 85% ANI:

| k (x) | y = 2 | y = 3 | y = 4 | y = 5 |
|---|---|---|---|---|
| 15 (72) | 247 | 165 | 110 | 74 |
| 16 (86) | 18.6 | 3.7 | 0.74 | 0.15 |
| 17 (102) | 1.4 | 0.08 | 0.005 | 3e-4 |
| 19 (142) | 0.008 | 4e-5 | 2e-7 | 1e-9 |

skani's k = 15 is safe only because its chaining requires coordinates on both genomes to agree.
Under the pure 1-D rule, k = 15 works only for ANI of about 95% and above. Two options:

- **(a) Keep the 1-D rule and raise k to 17–19.** skani stores seeds as `u32`, which allows
  k ≤ 16, so seeds would become `u64`. At 90% ANI and above, the variance cost is small: the
  SD of the ANI estimate barely depends on k there. At 85% ANI, though, the mean gap between
  matches grows from 1.4 kb to 2.0–2.7 kb, so boundaries are resolved less finely.
- **(b) Keep k = 15 and carry B's position with each match ("diagonal mode").** A match extends
  the current region only if its diagonal is within a band of the region's running diagonal.
  The diagonal is posA − posB, or posA + posB for reverse-strand matches. Any other match counts
  as a miss. This is still one linear scan, which makes it "chaining-lite". It removes spurious
  and paralogous matches, and it gives every region an interval on B as well (§8).

**Repeats.** Under the 1-D rule, every copy of a repeat in A matches if any copy exists in B.
This inflates AF for the genome with more repeats; IS elements and rRNA operons are typically
1–3% of a bacterial genome. Option (b) plus a 1-to-1 overlap filter removes the inflation.

## 7. Min-anchors (y)

y does two jobs.

**It rejects spurious regions.** Per §6, y = 3 gives fewer than 0.1 false regions per genome
at k ≥ 17, even with x ≈ 100.

**It sets the smallest usable contig.** At k = 17, the smallest contig expected to yield y
anchors (y·c/q) is:

| ANI | y = 3 | y = 4 | y = 5 |
|---|---|---|---|
| 0.99 | 0.4 kb | 0.6 kb | 0.7 kb |
| 0.95 | 0.9 kb | 1.2 kb | 1.5 kb |
| 0.90 | 2.2 kb | 3.0 kb | 3.7 kb |
| 0.85 | 5.9 kb | 7.9 kb | 9.9 kb |
| 0.80 | 16.7 kb | 22.2 kb | 27.8 kb |

y = 3 matches skani's minimum chain size. At 85% ANI, contigs shorter than about 6 kb add
nothing to AF whatever y is. This is inherent to c = 125; fragmented MAGs need a lower c.

## 8. The two-pass AF problem

**What skani does.** Its chains are approximately 1-to-1. It counts aligned bases on the query
once and divides that count by both genome lengths: `AF_ref = aligned query bases / |ref|`.
This is cheap, but correct only to the extent the chains really are 1-to-1.

**What we do.** The expensive step, finding shared k-mers, is symmetric. We intersect the two
sketches once:

1. Merge-join the two hash-sorted arrays. This is cache-friendly; skani's `get_anchors` also
   uses a galloping merge-join rather than hash lookups.
2. For each shared hash, set a bit for each of its entries in both genomes' position-ordered
   arrays, using the stored ranks.
3. Segment each genome with one linear scan over its bit vector.

The total cost is O(n_A + n_B). No 1-to-1 assumption is needed, and AF_A and AF_B are each
exact under the method's own definition.

In diagonal mode, each region already carries an interval on both genomes. A single scan then
gives both AFs, and overlapping projections onto B provide a 1-to-1 filter.

**Bonus: symmetric ANI.** Pooling gaps from both sides makes ANI(A, B) = ANI(B, A) exactly.
skani is asymmetric: it chooses which genome is the query using a score of length × contig
length.

## 9. Caveats shared with every k-mer ANI method

- **Heterogeneous divergence.** By Jensen's inequality, mean(a^k) > (mean a)^k. Pooled k-mer
  ANI therefore overestimates alignment ANI when divergence varies along the genome.
  Per-region estimates help. skani corrects this with a GBDT model trained against ANIm.
- **Indels** destroy k-mers, but alignment ANI ignores them, so k-mer ANI runs slightly low.
- **Fragmented MAGs.** Regions cannot cross contig ends (§7).
- **ANI floor** of about 80–85% at c = 125 (§3).

## 10. Proposed defaults (pending decisions in ROADMAP.md)

| Parameter | Flag | Default | Basis |
|---|---|---|---|
| k | `-k` | 17 (1-D) or 15 (diagonal mode) | §6 |
| c | `-c` | 125 | skani default, §3 |
| patience x | `--patience` | adaptive; integer override | §4 |
| per-gap false-break probability α | `--break-alpha` | 1e-3 | §4 |
| ANI floor for stage-1 patience | `--min-ani` | 0.85 (x ≈ 102 at k = 17, 72 at k = 15) | §3, §4 |
| min-anchors y | `--min-anchors` | 3 | §7 |
