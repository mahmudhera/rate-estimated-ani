use std::process::ExitCode;

use clap::Parser;
use kadence::cli::Cli;

fn main() -> ExitCode {
    // Usage errors exit with code 2 (clap's convention); runtime errors with code 1.
    let cli = Cli::parse();
    kadence::init_logging(cli.global.verbose);
    match kadence::run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err:#}");
            ExitCode::FAILURE
        }
    }
}
