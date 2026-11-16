# Roadmap: kadence

**kadence** estimates ANI and aligned fraction (AF) from positional FracMinHash sketches. It
finds homologous regions as runs of matched sketched k-mers, then estimates ANI from the
rate at which matches arrive inside those regions. It follows skani closely and replaces
skani's banded chaining with a run segmentation that has a statistical model behind it. See
[docs/THEORY.md](docs/THEORY.md) for the analysis behind every default below.

## Decisions

Open decisions (D1, D3 and D4 are needed before M3; D5 and D7 can wait):

| ID | Decision | Proposed default |
|---|---|---|
| D1 | Region finding: pure 1-D, or diagonal mode | Implement both, since they share about 90% of the code. Default chosen by M4 results. |
| D3 | Target ANI range | Accurate at ≥ 90%; reported down to an 85% floor |
| D4 | Patience | Adaptive (α = 1e-3, floor from `--min-ani 0.85`), with an integer `--patience` override |
| D5 | v1 CLI scope | `sketch` and `dist` only. `triangle` and `search` move to M5. |
| D7 | Main selling point against skani (speed, AF accuracy, or statistical rigor) | Decides benchmark priorities in M4 and M5 |

Decided (2026-11-16):

| ID | Decision | Outcome |
|---|---|---|
| D6 | Project and crate name, license | `kadence`; dual MIT OR Apache-2.0, like skani |
| D8 | CLI conventions | Our own minimal design (`-i`, `-o`, `-q`, `-r`, `--input-list`, `--query-list`, `--ref-list`), not a skani drop-in. The TSV output also uses our own schema. |
| D9 | Release artifacts | Static binaries on GitHub Releases for Linux (x86_64, aarch64, musl) and macOS (arm64, x86_64). crates.io and bioconda later. |
| D2 | Seed value | 64-bit `mm_hash64` hash of the canonical k-mer; k ≤ 31 (no 16-base limit as in skani). |
| D10 | Sketch file format | Custom versioned little-endian binary (`.ksk`) with a word-wise FNV-style checksum. bincode is unmaintained (RUSTSEC-2025-0141). |
| D11 | Sketch output layout | One file per genome, `OUT/<input file name>.ksk`. A consolidated database waits for M5. |
| D12 | Repeated k-mers | Keep every occurrence in the sketch; masking is a comparison-time option (M3). |
| D13 | Non-ACGT bases | Any non-ACGT base breaks k-mers (skani only breaks at N and reads other IUPAC codes as A). |
| D14 | Short contigs | Skip contigs < 500 bp (as skani), configurable with `--min-contig-len`. |
| D15 | SIMD | Port skani's AVX2 path only if the scalar path is slower than skani. It was, so it is ported, with runtime detection and a scalar fallback. |
| D16 | Output collisions | Fail before writing if two inputs map to the same sketch file; overwrite existing sketches. |

## Repository layout

Entries marked (M1) exist; the rest are planned.

```
Cargo.toml                 # one crate: lib + bin (M1)
src/main.rs, lib.rs        # entry point, run(), logging (M1)
src/cli.rs                 # clap derive definitions (M1)
src/params.rs              # SketchParams, MethodParams, provisional defaults (M1)
src/commands/{sketch,dist}.rs  # subcommands; stubs until M2/M3 (M1)
src/io.rs                  # FASTA/FASTQ via needletail, path lists (M2)
src/sketch/mod.rs          # Sketch (scan + set views), SketchBuilder, BitVec (M2)
src/sketch/kmer.rs         # rolling canonical k-mers, mm_hash64, scalar sampler (M2)
src/sketch/simd.rs         # AVX2 sampler, runtime-dispatched (M2)
src/sketch/format.rs       # .ksk binary format v1 (M2)
src/compare/{intersect,segment,estimate,af}.rs
src/output.rs              # TSV output (own schema)
tests/                     # CLI tests (M1); integration + simulation tests later
benches/                   # criterion micro-benchmarks
bench/                     # M4 Snakemake pipeline (skani, FastANI, MUMmer)
analysis/theory.py         # numbers in docs/THEORY.md
NOTICE                     # skani attribution (M2)
.github/workflows/{ci,release}.yml, .github/dependabot.yml   (M1)
deny.toml, rustfmt.toml, LICENSE-MIT, LICENSE-APACHE        (M1)
```

---

## M1: CLI and CI/CD

**Goal:** an empty but complete tool skeleton, with green CI and downloadable binaries.

- [x] Cargo crate with lib and bin targets: edition 2024, MSRV 1.85, `rustfmt.toml`, lint
  levels in `[lints]`, and a `profiling` build profile (release plus debug symbols).
- [x] CLI (clap derive):
  - Subcommands: `sketch`, `dist`.
  - Global flags: `-t/--threads` (default 3), `-v/--verbose` (repeatable).
  - `sketch`: `-i/--input`, `--input-list`, `-o/--output DIR`.
  - `dist`: `-q/--query`, `--query-list`, `-r/--ref`, `--ref-list`, `-o/--output`
    (default stdout).
  - Sketch flags (both subcommands): `-k/--kmer-size` (1–31), `-c/--scale`.
  - Method flags: `--mode {1d,diagonal}`, `--patience <auto|N>`, `--min-anchors`,
    `--break-alpha`, `--min-ani`.
  - Exit codes: an invalid command line gives 2. Valid arguments reach a stub that reports
    "not implemented yet" and exits with 1.
- [x] CI (GitHub Actions, on PRs and pushes to `main`):
  - `cargo fmt --check`, `clippy -D warnings`.
  - `cargo test` on ubuntu-latest and macos-latest.
  - An MSRV (1.85) check, rustdoc with `-D warnings`, and `cargo-deny` (licenses,
    advisories, sources).
  - Caching via `Swatinem/rust-cache`.
- [x] CD (on tags `vX.Y.Z` matching `Cargo.toml`): `gh release create` with generated notes,
  then `taiki-e/upload-rust-binary-action` builds static binaries for x86_64 and aarch64 Linux
  (musl, via `cross`) and for macOS arm64 and x86_64. Each archive gets a `.sha256` file.
- [x] Dependabot for cargo and actions. README, LICENSE-MIT, LICENSE-APACHE.

**Done when:**

- [ ] CI is green on a PR. Verified locally so far: fmt, clippy, 18 tests, MSRV 1.85
  build, docs, and cargo-deny.
- [ ] Tag `v0.0.1` produces downloadable binaries.
- [x] `--help` documents every parameter, with the units of patience stated (sketched
  k-mers). Covered by `tests/cli.rs`.

## M2: skani's parsing, plus our sketch representation

**Goal:** `kadence sketch` writes sketches that hold both the position-ordered "scan" view and
the hash-ordered "set" view.

- [x] Adapted from skani (credited in file headers and `NOTICE`):
  - FASTA/FASTQ reading via needletail 0.7: plain, gzip, bzip2, xz or zstd. Every record is a
    contig, and contigs under `--min-contig-len` (default 500) are skipped.
  - Rolling canonical k-mers and `mm_hash64`, with FracMinHash keeping hash < 2^64/c.
  - The AVX2 sampling path, ported to 64-bit hashes and k ≤ 31, with runtime detection and
    a scalar fallback. `KADENCE_NO_SIMD=1` forces the scalar path.
  - rayon parallelism across files, and list-file input (`--input-list`).
- [x] Dropped for now: AAI/ORF code, marker sketches (they return for screening in M5), and
  the GBDT model.
- [x] Sketch layout (`src/sketch/mod.rs`), struct-of-arrays:
  ```
  contig_names, contig_lengths, contig_offsets   // offsets index into the scan view
  hashes: Vec<u64>, positions: Vec<u32>, forward: BitVec   // scan view, position order
  by_hash: Vec<u32>          // set view: entry indices sorted by (hash, index)
  ```
  This takes about 16 bytes per entry, 600 KB for a 4.7 Mbp E. coli at c = 125. Positions are
  0-based k-mer starts within their contig. `forward` records whether the forward k-mer is
  the canonical one, which diagonal mode needs. The comparison between `by_hash` and a
  separate contiguous `sorted_hash` array moves to M3, where the merge-join lives.
- [x] Repeats are kept, not masked (decision D12). Every occurrence is an entry, and copies
  form contiguous runs in `by_hash`.
- [x] Versioned binary format `.ksk` v1 (`src/sketch/format.rs`, decision D10):
  - Header, arrays and a word-wise FNV-style checksum.
  - Written to a temporary file and renamed into place, with no fsync. fsync made sketching
    1.4× slower and is very slow on network filesystems.
  - Readers reject bad magic, other versions, bad checksums, truncation and trailing bytes,
    and validate every invariant.
- [x] Tests: 43 unit and 12 end-to-end (55), all passing. The AVX2 tests run only on
  x86_64 CPUs with AVX2, such as the Linux CI runners.
  - The rolling sampler equals a brute-force reference for k from 1 to 31 and several c, with
    lowercase, N and IUPAC bases. The dispatched (AVX2) path equals the scalar path, checked
    inside every sampler test.
  - Reverse complement gives the same hashes, mirrored positions and flipped strands.
  - Positions are strictly increasing within each contig, and offsets, the set-view
    permutation and its ordering are consistent (`Sketch::validate`).
  - The sampling rate is within 5 sd of 1/c on 2 Mbp of random sequence.
  - Serialization round-trips, corruption is detected, and gzip input matches plain input.
  - CLI: list input, name collisions (nothing written), partial failures (good inputs still
    written, exit code 1), overwriting, empty input, `--min-contig-len`, and byte-identical
    output for 1 and 4 threads.
  - **Changed: no seed-set parity test with skani.** skani v0.4.0's seed is not the canonical
    k-mer at its position. Both its scalar and AVX2 paths take the canonical form as
    min(forward of the last k-mer, reverse complement of the first k-mer) of a 21-base
    window. Ported to Python, 66.6% of skani's seeds differ from the canonical 15-mer at
    their position. The seeds are still strand-symmetric and survive mutation like k ≈ 15–16
    (measured 15.1–15.9), but an exact seed-set match is impossible by design. The
    brute-force reference replaces the parity test. Worth reporting upstream.

**Benchmark** (skani's 7 bacterial test genomes × 4 = 28 files, ~140 Mbp; Xeon Gold 5416S
shared with other jobs at load ≈ 30/64, so medians of interleaved runs):

Medians of 7–9 interleaved skani/kadence pairs, at k = 15 and c = 125 to match skani.
Sessions A and B were about 10 minutes apart.

| Run | skani | kadence | kadence / skani |
|---|---|---|---|
| 1 thread, scalar sampling only | 563 ms CPU | 713 ms CPU | 1.27 CPU |
| 1 thread, AVX2 (session A) | 586 ms CPU, 519 ms wall | 586 ms CPU, 587 ms wall | 1.00 CPU, 1.14 wall |
| 1 thread, final build (session B, load ≈ 33) | 841 ms CPU, 778 ms wall | 1081 ms CPU, 1081 ms wall | 1.22 CPU, 1.33 wall (pairs ranged 0.86–1.72) |
| 8 threads, AVX2 (session A) | 151 ms wall | 171 ms wall | 1.18 wall |
| 8 threads, final build (session B) | 191 ms wall, 1117 ms CPU | 195 ms wall, 1206 ms CPU | 1.05 wall, 1.08 CPU |

Measured within kadence:

- Dropping fsync made it 1.4× faster.
- AVX2 sampling is 17% faster end-to-end than scalar sampling, and about 2× faster for the
  sampling stage alone.
- Two other ideas were tried and dropped because they were slower: a two-stream scalar loop
  (2.30 vs 2.22 ns/base) and three-pass blocking (about 1.5× slower).

Stage costs (scalar, quieter moment), in ns per base:

| Stage | ns/base |
|---|---|
| Parse + sequence copy | 1.1 |
| Sampling | 2.2 |
| Entry push + set-view sort | 0.4 |
| Serialization | 0.1 |

**Verdict:** within about 5–20% of skani, and not reproducibly faster. The machine was shared
(another user's job used about 29 cores), and identical runs varied by up to 1.7×. The user
time at 8 threads (734 vs 648 ms) suggests kadence still does about 13% more per-core work.
Next steps:

- Re-measure on an idle node with hyperfine, as part of M4.
- If the gap holds, candidates for M5: stream multi-line FASTA without copying the sequence,
  pre-size the entry vectors, and overlap writing with sketching.

**Done when:**

- [x] Sketches build for skani's E. coli fixtures.
- [x] All invariants pass.
- [ ] Sketching throughput is at least skani's. Near parity, but not shown conclusively (see
  the benchmark above); re-measure on an idle node in M4.

## M3: method v1 and testing

**Goal:** `kadence dist` reports ANI, AF_A, AF_B and a confidence interval, from FASTA or
sketch input.

- [ ] `intersect`: galloping merge-join over the hash-sorted arrays. It produces
  `matched_A` and `matched_B` bit vectors, and in diagonal mode `(rankA, posB, strand)` anchors.
- [ ] `segment`: one linear scan per contig over the bit vectors, applying patience and
  min-anchors. Contig ends close a region. Diagonal mode adds a band test against the
  region's running diagonal.
- [ ] `estimate`:
  - Truncated-geometric MLE of q over the pooled inter-match gaps of both genomes, solved
    by bisection.
  - ANI = q̂^(1/k), with a Fisher-information standard error and a delta-method CI.
  - The naive m/n estimate is also kept, for diagnostics.
- [ ] Adaptive patience loop: start from floor x, estimate q̂, set x(q̂, α), and repeat until
  x is stable (at most 3 iterations).
- [ ] AF: the summed base-pair span of accepted regions, plus boundary padding of the expected
  half-gap (capped at contig ends), divided by genome length. Do this for each genome. Add a
  flag that disables padding, for ablation.
- [ ] Output: our own documented TSV schema (decision D8). Columns: query and reference
  names, ANI, `ANI_lo`, `ANI_hi`, AF of the query, AF of the reference, `q_hat`,
  `patience_used`, `n_regions`, and `n_anchors`. The M4 harness maps these columns to
  skani's.
- [ ] Tests:
  - Unit tests:
    - Segmenter edge cases: empty input, all matches, x = 1, y = 1, contig boundaries,
      and strand.
    - MLE round-trip on simulated gaps.
  - Property tests (proptest):
    - ANI(A,B) = ANI(B,A).
    - 0 ≤ AF ≤ 1.
    - Self-comparison gives ANI = AF = 1.
    - Results are unchanged when contigs are reordered or reverse-complemented.
  - Simulation tests, run in CI at small scale:
    - Build a random genome, then mutate it to target ANI from 0.80 to 0.999, with indels.
    - Insert islands of known size and fragment the result to a target N50.
    - Check ANI and AF against the truth, within tolerances set per ANI bin.
  - Golden tests: skani's E. coli fixtures, compared approximately with skani's output as a
    sanity check.

**Done when:**

- All tests are green in CI.
- The simulation error is within tolerance for ANI ≥ 85%.
- The output schema is documented.

## M4: benchmark against skani

**Goal:** a reproducible report showing where kadence matches, beats or trails skani.

- [ ] Comparators:
  - skani: default (c = 125) and `--slow` (c = 30).
  - FastANI.
  - MUMmer4 `dnadiff` (ANIm), as ground truth for both ANI and AF.
  - Optionally pyani ANIb.
- [ ] Datasets:
  - (1) A simulated grid over ANI, AF, N50 and genome size.
  - (2) Real bacterial pairs spanning 80–100% ANI.
  - (3) MAG sets (fragmented, incomplete), following the skani paper.
  - (4) The EvANI benchmark.
- [ ] Metrics:
  - ANI bias and MAE against ANIm, per ANI bin, plus Pearson and Spearman correlation.
  - AF MAE against dnadiff.
  - Precision and recall for the 95% species boundary.
  - Runtime and peak RSS for sketching and for all-vs-all `dist` on about 1k genomes.
- [ ] Ablations: 1-D vs diagonal mode, fixed vs adaptive patience, naive vs MLE estimator,
  and k ∈ {15, 17, 19}.
- [ ] Infrastructure: a Snakemake pipeline in `bench/`, with a conda env for the comparators,
  producing plots and a report.

**Done when:** the report exists, and D1, D4 and the default k are settled by data.

## M5: improve, reimplement, retest (iterative)

Candidates, ranked by expected impact. M4 results will reorder them.

1. **Two-state HMM segmentation.** It replaces the patience and min-anchors heuristic with
   Viterbi or posterior decoding, adapts to local ANI, and gives boundary uncertainty
   (THEORY §4).
2. **Diagonal mode with a 1-to-1 overlap filter**, if M4 shows the 1-D mode suffers from
   repeats or spurious matches.
3. **Aggregation choice:** pooled, length-weighted per-region, median, or trimmed. Pick the
   one that minimizes Jensen bias against ANIm, and add a small calibration only if needed.
4. **Indel correction**, for k-mers lost to indels.
5. **Scale-out:** `triangle` and `search`, with marker-based screening, an inverted index,
   and a consolidated database (skani's design).
6. **Performance:** SIMD merge-join (skani has AVX2), mmap zero-copy sketches, and better
   rayon scheduling over pairs.
7. **Lower c, chosen automatically**, for small genomes and fragmented MAGs.

Each iteration: change, retest (M3 suite), rebenchmark (M4 pipeline), and record the result
in `bench/RESULTS.md`.
