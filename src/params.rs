//! Parameter types shared by the command line and the library.
//!
//! The method defaults are provisional. They come from `docs/THEORY.md` (section 10)
//! and will be revisited once the open decisions in `ROADMAP.md` are settled by benchmarks.

use std::fmt;
use std::num::NonZeroU32;
use std::str::FromStr;

/// Default k-mer size (provisional; THEORY section 6).
pub const DEFAULT_K: u8 = 17;
/// Largest k whose 2-bit encoding fits in a `u64`.
pub const MAX_K: u8 = 31;
/// Default FracMinHash scale: keep about one in `c` k-mers.
pub const DEFAULT_C: u32 = 125;
/// Contigs shorter than this are left out of a sketch and of the genome length (as in skani).
pub const DEFAULT_MIN_CONTIG_LEN: u32 = 500;
/// Default number of matched sketched k-mers required to accept a region.
pub const DEFAULT_MIN_ANCHORS: u32 = 3;
/// Default per-gap false-break probability for adaptive patience.
pub const DEFAULT_BREAK_ALPHA: f64 = 1e-3;
/// Default lowest ANI of interest, which sets the starting patience.
pub const DEFAULT_MIN_ANI: f64 = 0.85;
/// Default number of worker threads.
pub const DEFAULT_THREADS: usize = 3;

/// How many consecutive non-matching sketched k-mers a region tolerates before it closes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Patience {
    /// Derive patience from the pair's own match rate (THEORY section 4).
    #[default]
    Auto,
    /// Close a region after this many consecutive non-matching sketched k-mers.
    Fixed(NonZeroU32),
}

impl fmt::Display for Patience {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Auto => f.write_str("auto"),
            Self::Fixed(n) => write!(f, "{n}"),
        }
    }
}

impl FromStr for Patience {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.eq_ignore_ascii_case("auto") {
            return Ok(Self::Auto);
        }
        let n: u32 = s
            .parse()
            .map_err(|_| format!("expected `auto` or a positive integer, got `{s}`"))?;
        NonZeroU32::new(n)
            .map(Self::Fixed)
            .ok_or_else(|| "patience must be at least 1".to_string())
    }
}

/// How matched k-mers are grouped into homologous regions.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Mode {
    /// Runs of matches along one genome's sketch.
    #[default]
    #[value(name = "1d")]
    OneD,
    /// Runs of matches that also stay on a consistent diagonal in the other genome.
    Diagonal,
}

impl fmt::Display for Mode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::OneD => "1d",
            Self::Diagonal => "diagonal",
        })
    }
}

/// Parameters that define a sketch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SketchParams {
    /// k-mer size, `1..=MAX_K`.
    pub k: u8,
    /// FracMinHash scale; about one in `c` k-mers is kept.
    pub c: u32,
    /// Contigs shorter than this many bases are skipped.
    pub min_contig_len: u32,
}

impl SketchParams {
    /// Whether two sketches sample the same k-mers, so their hashes can be compared.
    /// `min_contig_len` does not affect which k-mers are sampled.
    pub fn compatible_with(&self, other: &Self) -> bool {
        self.k == other.k && self.c == other.c
    }
}

impl Default for SketchParams {
    fn default() -> Self {
        Self {
            k: DEFAULT_K,
            c: DEFAULT_C,
            min_contig_len: DEFAULT_MIN_CONTIG_LEN,
        }
    }
}

/// Parameters for finding homologous regions and estimating ANI and AF.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MethodParams {
    /// How matches are grouped into regions.
    pub mode: Mode,
    /// Consecutive misses tolerated inside a region.
    pub patience: Patience,
    /// Matches required to accept a region.
    pub min_anchors: NonZeroU32,
    /// Per-gap false-break probability, in `(0, 1)`.
    pub break_alpha: f64,
    /// Lowest ANI of interest, in `(0, 1)`.
    pub min_ani: f64,
}

impl Default for MethodParams {
    fn default() -> Self {
        Self {
            mode: Mode::default(),
            patience: Patience::default(),
            min_anchors: NonZeroU32::new(DEFAULT_MIN_ANCHORS).expect("default is non-zero"),
            break_alpha: DEFAULT_BREAK_ALPHA,
            min_ani: DEFAULT_MIN_ANI,
        }
    }
}

/// Parses a probability strictly between 0 and 1.
pub fn parse_open_unit_interval(s: &str) -> Result<f64, String> {
    let x: f64 = s
        .parse()
        .map_err(|_| format!("expected a number, got `{s}`"))?;
    if x > 0.0 && x < 1.0 {
        Ok(x)
    } else {
        Err(format!("must be strictly between 0 and 1, got {x}"))
    }
}

/// Parses a positive integer (at least 1).
pub fn parse_positive_usize(s: &str) -> Result<usize, String> {
    match s.parse::<usize>() {
        Ok(0) => Err("must be at least 1".to_string()),
        Ok(n) => Ok(n),
        Err(_) => Err(format!("expected a positive integer, got `{s}`")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn patience_parses_auto_and_integers() {
        assert_eq!("auto".parse::<Patience>(), Ok(Patience::Auto));
        assert_eq!("AUTO".parse::<Patience>(), Ok(Patience::Auto));
        assert_eq!(
            "12".parse::<Patience>(),
            Ok(Patience::Fixed(NonZeroU32::new(12).unwrap()))
        );
    }

    #[test]
    fn patience_rejects_zero_and_garbage() {
        assert!("0".parse::<Patience>().is_err());
        assert!("-3".parse::<Patience>().is_err());
        assert!("soon".parse::<Patience>().is_err());
    }

    #[test]
    fn patience_display_round_trips() {
        for p in [Patience::Auto, Patience::Fixed(NonZeroU32::new(7).unwrap())] {
            assert_eq!(p.to_string().parse::<Patience>(), Ok(p));
        }
    }

    #[test]
    fn open_unit_interval_bounds() {
        assert_eq!(parse_open_unit_interval("0.85"), Ok(0.85));
        assert_eq!(parse_open_unit_interval("1e-3"), Ok(1e-3));
        for bad in ["0", "1", "1.5", "-0.1", "NaN", "x"] {
            assert!(parse_open_unit_interval(bad).is_err(), "{bad} accepted");
        }
    }

    #[test]
    fn positive_usize_bounds() {
        assert_eq!(parse_positive_usize("8"), Ok(8));
        assert!(parse_positive_usize("0").is_err());
        assert!(parse_positive_usize("-1").is_err());
    }

    #[test]
    fn defaults_are_valid() {
        let m = MethodParams::default();
        assert!(m.break_alpha > 0.0 && m.break_alpha < 1.0);
        assert!(m.min_ani > 0.0 && m.min_ani < 1.0);
        let s = SketchParams::default();
        assert!((1..=MAX_K).contains(&s.k));
        assert!(s.c >= 1);
    }
}
