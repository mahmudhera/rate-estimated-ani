# Roadmap: kadence

**kadence** estimates ANI and aligned fraction (AF) from positional FracMinHash sketches. It
finds homologous regions as runs of matched sketched k-mers, then estimates ANI from the
rate at which matches arrive inside those regions. It follows skani closely and replaces
skani's banded chaining with a run segmentation that has a statistical model behind it. See
[docs/THEORY.md](docs/THEORY.md) for the analysis behind every default below.

## Decisions

Open decisions (D1–D4 are needed before M3; D5 and D7 can wait):

| ID | Decision | Proposed default |
|---|---|---|
| D1 | Region finding: pure 1-D, or diagonal mode | Implement both, since they share about 90% of the code. Default chosen by M4 results. |
| D2 | Seed type | `u64` hashes (k ≤ 31). This drops skani's k ≤ 16 limit. |
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

## Repository layout

Entries marked (M1) exist; the rest are planned.

```
Cargo.toml                 # one crate: lib + bin (M1)
src/main.rs, lib.rs        # entry point, run(), logging (M1)
src/cli.rs                 # clap derive definitions (M1)
src/params.rs              # SketchParams, MethodParams, provisional defaults (M1)
src/commands/{sketch,dist}.rs  # subcommands; stubs until M2/M3 (M1)
src/io/fastx.rs            # adapted from skani file_io.rs (needletail)
src/sketch/{seeding,hash,format,mod}.rs
src/compare/{intersect,segment,estimate,af}.rs
src/output.rs              # TSV output (own schema)
tests/                     # CLI tests (M1); integration + simulation tests later
benches/                   # criterion micro-benchmarks
bench/                     # M4 Snakemake pipeline (skani, FastANI, MUMmer)
analysis/theory.py         # numbers in docs/THEORY.md
NOTICE                     # skani attribution (added in M2, with the adapted code)
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

- [ ] Adapt from skani, crediting it in file headers and in a top-level `NOTICE`:
  - FASTA/FASTQ(.gz) reading via needletail, with minimum contig length 500 and N handling.
  - Rolling canonical k-mers and `mm_hash64`.
  - The AVX2 seeding fast path, ported to `u64`.
  - rayon parallelism across files, and list-file input.
- [ ] Drop for now: AAI/ORF code, marker sketches (they return for screening in M5), and the
  GBDT model.
- [ ] Sketch layout, struct-of-arrays:
  ```
  contig_names, contig_lens, contig_offsets   // offsets index into the pos-ordered arrays
  hash: Vec<u64>, pos: Vec<u32>, strand: bitvec    // position order (scan view)
  by_hash: Vec<u32>                                 // ranks sorted by (hash, rank) (set view)
  ```
  This takes about 16 bytes per seed, roughly 650 KB per 5 Mbp genome at c = 125. Benchmark
  `by_hash` against a separate contiguous `sorted_hash: Vec<u64>` to see which merge-join is
  faster.
- [ ] Repeat masking: drop seeds that occur more than 2500/c times (skani's rule), configurable.
- [ ] Versioned binary format: start with bincode, as skani does, and evaluate mmap or rkyv
  zero-copy loading in M5.
- [ ] Tests:
  - Hand-checked k-mers on tiny sequences.
  - Reverse-complement invariance: same hash multiset, mirrored positions.
  - Positions monotone within each contig, with contig offsets consistent.
  - Statistical check that the sampling rate is 1/c on random sequence.
  - Serialization round-trip.
  - Seed-set parity with skani at k = 15, as a sanity check.

**Done when:**

- Sketches build for skani's E. coli fixtures.
- All invariants pass.
- Sketching throughput is at least skani's, measured with hyperfine.

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
