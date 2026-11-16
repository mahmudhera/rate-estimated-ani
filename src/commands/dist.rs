//! `kadence dist`.

use anyhow::{Result, bail};

use crate::cli::DistArgs;
use crate::params::{MethodParams, SketchParams};

/// Estimates ANI and AF between queries and references. Not implemented yet (milestone M3).
pub fn run(args: &DistArgs, threads: usize) -> Result<()> {
    let sketch = SketchParams::from(&args.sketch);
    let method = MethodParams::from(&args.method);
    log::debug!(
        "dist: {} query file(s), query list {:?}, {} reference file(s), reference list {:?}, \
         output {:?}, {sketch:?}, {method:?}, {threads} thread(s)",
        args.query.len(),
        args.query_list,
        args.reference.len(),
        args.reference_list,
        args.output,
    );
    bail!("`kadence dist` is not implemented yet (planned for milestone M3)")
}
