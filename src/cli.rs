//! Command-line interface definition.

use std::num::NonZeroU32;
use std::path::PathBuf;

use clap::{ArgGroup, Args, Parser, Subcommand};

use crate::params::{self, MethodParams, Mode, Patience, SketchParams};

/// Rate-based ANI and aligned fraction from positional FracMinHash sketches.
#[derive(Debug, Parser)]
#[command(name = "kadence", version, about, propagate_version = true)]
pub struct Cli {
    #[command(flatten)]
    pub global: GlobalArgs,

    #[command(subcommand)]
    pub command: Command,
}

/// Options accepted by every subcommand.
#[derive(Debug, Args)]
pub struct GlobalArgs {
    /// Number of worker threads.
    #[arg(
        short,
        long,
        global = true,
        help_heading = "Global options",
        default_value_t = params::DEFAULT_THREADS,
        value_parser = params::parse_positive_usize
    )]
    pub threads: usize,

    /// More log output (-v debug, -vv trace). RUST_LOG overrides this.
    #[arg(
        short,
        long,
        global = true,
        help_heading = "Global options",
        action = clap::ArgAction::Count
    )]
    pub verbose: u8,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Sketch genomes into FracMinHash sketches that keep k-mer positions.
    Sketch(SketchArgs),
    /// Estimate ANI and aligned fraction between query and reference genomes.
    Dist(DistArgs),
}

#[derive(Debug, Args)]
#[command(group(
    ArgGroup::new("inputs")
        .required(true)
        .multiple(true)
        .args(["input", "input_list"])
))]
pub struct SketchArgs {
    /// FASTA files to sketch (gzip allowed).
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    pub input: Vec<PathBuf>,

    /// File listing FASTA paths, one per line.
    #[arg(long, value_name = "FILE")]
    pub input_list: Option<PathBuf>,

    /// Directory to write sketches into.
    #[arg(short, long, value_name = "DIR")]
    pub output: PathBuf,

    #[command(flatten)]
    pub sketch: SketchOpts,
}

#[derive(Debug, Args)]
#[command(group(
    ArgGroup::new("queries")
        .required(true)
        .multiple(true)
        .args(["query", "query_list"])
))]
#[command(group(
    ArgGroup::new("references")
        .required(true)
        .multiple(true)
        .args(["reference", "reference_list"])
))]
pub struct DistArgs {
    /// Query genomes: FASTA files (gzip allowed) or kadence sketches.
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    pub query: Vec<PathBuf>,

    /// File listing query paths, one per line.
    #[arg(long, value_name = "FILE")]
    pub query_list: Option<PathBuf>,

    /// Reference genomes: FASTA files (gzip allowed) or kadence sketches.
    #[arg(short, long = "ref", num_args = 1.., value_name = "FILE")]
    pub reference: Vec<PathBuf>,

    /// File listing reference paths, one per line.
    #[arg(long = "ref-list", value_name = "FILE")]
    pub reference_list: Option<PathBuf>,

    /// Output TSV file [default: stdout].
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Used when an input is FASTA rather than a sketch.
    #[command(flatten)]
    pub sketch: SketchOpts,

    #[command(flatten)]
    pub method: MethodOpts,
}

/// Sketch parameters.
#[derive(Debug, Args)]
#[command(next_help_heading = "Sketch options")]
pub struct SketchOpts {
    /// k-mer size (1-31).
    #[arg(
        short,
        long = "kmer-size",
        value_name = "K",
        default_value_t = params::DEFAULT_K,
        value_parser = clap::value_parser!(u8).range(1..=i64::from(params::MAX_K))
    )]
    pub k: u8,

    /// FracMinHash scale: keep about one in C k-mers.
    #[arg(
        short,
        long = "scale",
        value_name = "C",
        default_value_t = params::DEFAULT_C,
        value_parser = clap::value_parser!(u32).range(1..)
    )]
    pub c: u32,

    /// Skip contigs shorter than N bases; they also don't count toward genome length.
    #[arg(long, value_name = "N", default_value_t = params::DEFAULT_MIN_CONTIG_LEN)]
    pub min_contig_len: u32,
}

/// Region-finding and estimation parameters. Defaults are provisional (docs/THEORY.md).
#[derive(Debug, Args)]
#[command(next_help_heading = "Method options (provisional defaults, see docs/THEORY.md)")]
pub struct MethodOpts {
    /// How matched k-mers are grouped into homologous regions.
    #[arg(long, value_enum, default_value_t = Mode::OneD)]
    pub mode: Mode,

    /// Consecutive non-matching sketched k-mers before a region closes, or `auto` to
    /// derive it from the pair's match rate.
    #[arg(long, value_name = "auto|N", default_value_t = Patience::Auto)]
    pub patience: Patience,

    /// Matched sketched k-mers required to accept a region.
    #[arg(
        long,
        value_name = "N",
        default_value_t = params::DEFAULT_MIN_ANCHORS,
        value_parser = clap::value_parser!(u32).range(1..)
    )]
    pub min_anchors: u32,

    /// Per-gap false-break probability used by `--patience auto`.
    #[arg(
        long,
        value_name = "P",
        default_value_t = params::DEFAULT_BREAK_ALPHA,
        value_parser = params::parse_open_unit_interval
    )]
    pub break_alpha: f64,

    /// Lowest ANI of interest (0-1); sets the starting patience for `--patience auto`.
    #[arg(
        long,
        value_name = "ANI",
        default_value_t = params::DEFAULT_MIN_ANI,
        value_parser = params::parse_open_unit_interval
    )]
    pub min_ani: f64,
}

impl From<&SketchOpts> for SketchParams {
    fn from(o: &SketchOpts) -> Self {
        Self {
            k: o.k,
            c: o.c,
            min_contig_len: o.min_contig_len,
        }
    }
}

impl From<&MethodOpts> for MethodParams {
    fn from(o: &MethodOpts) -> Self {
        Self {
            mode: o.mode,
            patience: o.patience,
            // The value parser guarantees min_anchors >= 1.
            min_anchors: NonZeroU32::new(o.min_anchors).expect("validated by clap"),
            break_alpha: o.break_alpha,
            min_ani: o.min_ani,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;
    use clap::error::ErrorKind;

    fn parse(args: &[&str]) -> Result<Cli, clap::Error> {
        Cli::try_parse_from(std::iter::once("kadence").chain(args.iter().copied()))
    }

    #[test]
    fn cli_definition_is_consistent() {
        Cli::command().debug_assert();
    }

    #[test]
    fn sketch_requires_input_and_output() {
        let e = parse(&["sketch", "-o", "out"]).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::MissingRequiredArgument);
        let e = parse(&["sketch", "-i", "a.fa"]).unwrap_err();
        assert_eq!(e.kind(), ErrorKind::MissingRequiredArgument);
    }

    #[test]
    fn sketch_accepts_files_and_list() {
        let cli = parse(&[
            "sketch",
            "-i",
            "a.fa",
            "b.fa",
            "--input-list",
            "l.txt",
            "-o",
            "d",
        ])
        .unwrap();
        let Command::Sketch(args) = cli.command else {
            panic!("expected sketch")
        };
        assert_eq!(args.input.len(), 2);
        assert_eq!(args.input_list, Some(PathBuf::from("l.txt")));
        assert_eq!(SketchParams::from(&args.sketch), SketchParams::default());
    }

    #[test]
    fn dist_requires_both_sides() {
        for args in [
            &["dist", "-q", "a.fa"][..],
            &["dist", "-r", "b.fa"][..],
            &["dist", "--query-list", "q.txt"][..],
        ] {
            let e = parse(args).unwrap_err();
            assert_eq!(e.kind(), ErrorKind::MissingRequiredArgument, "{args:?}");
        }
    }

    #[test]
    fn dist_defaults_match_params_defaults() {
        let cli = parse(&["dist", "-q", "a.fa", "-r", "b.fa"]).unwrap();
        let Command::Dist(args) = cli.command else {
            panic!("expected dist")
        };
        assert_eq!(MethodParams::from(&args.method), MethodParams::default());
        assert_eq!(SketchParams::from(&args.sketch), SketchParams::default());
        assert_eq!(cli.global.threads, params::DEFAULT_THREADS);
        assert!(args.output.is_none());
    }

    #[test]
    fn dist_parses_method_options() {
        let cli = parse(&[
            "dist",
            "-q",
            "a.fa",
            "b.fa",
            "--ref-list",
            "r.txt",
            "--mode",
            "diagonal",
            "--patience",
            "40",
            "--min-anchors",
            "4",
            "--break-alpha",
            "0.01",
            "--min-ani",
            "0.9",
            "-k",
            "15",
            "-c",
            "200",
            "-t",
            "8",
        ])
        .unwrap();
        let Command::Dist(args) = cli.command else {
            panic!("expected dist")
        };
        let m = MethodParams::from(&args.method);
        assert_eq!(m.mode, Mode::Diagonal);
        assert_eq!(m.patience, Patience::Fixed(NonZeroU32::new(40).unwrap()));
        assert_eq!(m.min_anchors.get(), 4);
        assert_eq!(m.break_alpha, 0.01);
        assert_eq!(m.min_ani, 0.9);
        assert_eq!(
            SketchParams::from(&args.sketch),
            SketchParams {
                k: 15,
                c: 200,
                min_contig_len: params::DEFAULT_MIN_CONTIG_LEN
            }
        );
        assert_eq!(args.query.len(), 2);
        assert_eq!(cli.global.threads, 8);
    }

    #[test]
    fn out_of_range_values_are_rejected() {
        let base = ["dist", "-q", "a.fa", "-r", "b.fa"];
        for extra in [
            &["-k", "0"][..],
            &["-k", "32"][..],
            &["-c", "0"][..],
            &["--patience", "0"][..],
            &["--min-anchors", "0"][..],
            &["--break-alpha", "1"][..],
            &["--min-ani", "1.2"][..],
            &["--mode", "2d"][..],
            &["-t", "0"][..],
        ] {
            let args: Vec<&str> = base.iter().chain(extra).copied().collect();
            assert!(parse(&args).is_err(), "{extra:?} accepted");
        }
    }
}
