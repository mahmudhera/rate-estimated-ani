# kadence

Estimate average nucleotide identity (ANI) and aligned fraction (AF) between genomes from
positional FracMinHash sketches.

> **Status: pre-alpha.** Milestone M1 (CLI and CI/CD) is in place. `sketch` and `dist` validate
> their arguments but do not compute anything yet. See [ROADMAP.md](ROADMAP.md).

kadence follows [skani](https://github.com/bluenote-1577/skani) closely. Each sketch keeps the
positions of its sampled k-mers. When two genomes are compared, runs of matching sketched k-mers
mark homologous regions. ANI comes from the rate at which matches arrive inside those regions,
and AF from how much of each genome the regions cover. The theory behind the method and its
defaults is in [docs/THEORY.md](docs/THEORY.md) ([PDF](docs/THEORY.pdf)).

## Install

**Prebuilt binaries.** Each [GitHub release](../../releases) includes static Linux builds
(x86_64, aarch64) and macOS builds (Apple silicon, Intel), each with a `.sha256` checksum.

**From source** (Rust 1.85 or newer):

```sh
cargo install --path .
# or: cargo build --release   # binary at target/release/kadence
```

## Usage (planned interface)

```sh
kadence sketch -i genomes/*.fa -o sketches/
kadence sketch --input-list genomes.txt -o sketches/

kadence dist -q a.fa -r b.fa
kadence dist -q sketches/a.sketch -r refs/*.fa -o results.tsv
kadence dist --query-list queries.txt --ref-list refs.txt -t 16
```

`dist` accepts FASTA files (gzip allowed) or kadence sketches. Run `kadence <command> --help`
for every option.

Method options have **provisional** defaults; benchmarks in milestone M4 will settle them.

| Option | Meaning | Default |
|---|---|---|
| `--mode <1d\|diagonal>` | Group matches along one genome, or also require a consistent diagonal in the other | `1d` |
| `--patience <auto\|N>` | Consecutive non-matching sketched k-mers before a region closes; `auto` derives it from the pair's match rate | `auto` |
| `--min-anchors <N>` | Matched sketched k-mers required to accept a region | `3` |
| `--break-alpha <P>` | Per-gap false-break probability used by `--patience auto` | `0.001` |
| `--min-ani <ANI>` | Lowest ANI of interest; sets the starting patience | `0.85` |
| `-k, --kmer-size <K>` | k-mer size (1–31) | `17` |
| `-c, --scale <C>` | FracMinHash scale: keep about one in C k-mers | `125` |

Exit codes: `0` success, `1` runtime error, `2` invalid command line.

## Development

```sh
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo build --profile profiling   # release speed with debug symbols, for profiling
```

CI runs on every pull request and on pushes to `main`. It runs rustfmt, clippy, tests on Linux
and macOS, a check with the minimum supported Rust version (1.85), docs, and `cargo-deny` for
licenses and advisories.

**Releasing:** bump `version` in `Cargo.toml`, merge to `main`, then push a matching tag:

```sh
git tag v0.0.1 && git push origin v0.0.1
```

The release workflow builds the binaries and publishes them on the GitHub release.

## License

Licensed under either of [Apache License 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT),
at your option.

## Acknowledgements

kadence builds on skani (J. Shaw and Y. W. Yu, *Nature Methods* 2023). Starting with milestone
M2, it adapts skani's sequence parsing and sketching code, which is available under MIT OR
Apache-2.0.
