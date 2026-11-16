//! `kadence sketch`.

use anyhow::{Result, bail};

use crate::cli::SketchArgs;
use crate::params::SketchParams;

/// Sketches the input genomes. Not implemented yet (milestone M2).
pub fn run(args: &SketchArgs, threads: usize) -> Result<()> {
    let params = SketchParams::from(&args.sketch);
    log::debug!(
        "sketch: {} input file(s), input list {:?}, output {}, {params:?}, {threads} thread(s)",
        args.input.len(),
        args.input_list,
        args.output.display(),
    );
    bail!("`kadence sketch` is not implemented yet (planned for milestone M2)")
}
