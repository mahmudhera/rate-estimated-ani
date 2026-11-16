//! kadence: rate-based ANI and aligned fraction from positional FracMinHash sketches.
//!
//! See `docs/THEORY.md` for the method and `ROADMAP.md` for the implementation plan.

pub mod cli;
pub mod commands;
pub mod params;

use cli::{Cli, Command};

/// Runs the subcommand selected on the command line.
pub fn run(cli: &Cli) -> anyhow::Result<()> {
    match &cli.command {
        Command::Sketch(args) => commands::sketch::run(args, cli.global.threads),
        Command::Dist(args) => commands::dist::run(args, cli.global.threads),
    }
}

/// Sets up logging to stderr: info by default, `-v` for debug, `-vv` for trace.
/// `RUST_LOG`, when set, takes precedence.
pub fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => log::LevelFilter::Info,
        1 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    env_logger::Builder::new()
        .filter_level(level)
        .parse_default_env()
        .init();
}
