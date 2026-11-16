//! `kadence sketch`: sketch FASTA/FASTQ files into one `.ksk` file per genome.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::cli::SketchArgs;
use crate::io;
use crate::params::SketchParams;
use crate::sketch::format;

/// One input file and the sketch file it produces.
#[derive(Debug, PartialEq, Eq)]
struct Job {
    input: PathBuf,
    output: PathBuf,
}

/// What one sketched genome contributed to the summary.
struct Done {
    contigs: usize,
    bases: u64,
    kmers: usize,
}

/// Sketches every input in parallel and writes `OUT/<input file name>.ksk` for each.
///
/// All inputs are attempted. A failing input is reported, and the command fails at the end.
pub fn run(args: &SketchArgs, threads: usize) -> Result<()> {
    let params = SketchParams::from(&args.sketch);
    let jobs = plan_jobs(collect_inputs(args)?, &args.output)?;
    fs::create_dir_all(&args.output).with_context(|| {
        format!(
            "failed to create output directory {}",
            args.output.display()
        )
    })?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("failed to start worker threads")?;
    log::info!(
        "sketching {} file(s) with k = {}, c = {}, min contig length {} bp, {threads} thread(s)",
        jobs.len(),
        params.k,
        params.c,
        params.min_contig_len,
    );

    let start = Instant::now();
    let results: Vec<Result<Done>> =
        pool.install(|| jobs.par_iter().map(|job| sketch_one(job, params)).collect());

    let (mut genomes, mut contigs, mut bases, mut kmers, mut failures) = (0, 0, 0u64, 0, 0);
    for result in results {
        match result {
            Ok(d) => {
                genomes += 1;
                contigs += d.contigs;
                bases += d.bases;
                kmers += d.kmers;
            }
            Err(e) => {
                log::error!("{e:#}");
                failures += 1;
            }
        }
    }
    log::info!(
        "sketched {genomes} genome(s): {contigs} contigs, {bases} bp, {kmers} sampled k-mers \
         in {:.2} s; sketches are in {}",
        start.elapsed().as_secs_f64(),
        args.output.display(),
    );
    if failures > 0 {
        bail!("{failures} of {} input file(s) failed", jobs.len());
    }
    Ok(())
}

/// Inputs from `--input`, then from `--input-list`.
fn collect_inputs(args: &SketchArgs) -> Result<Vec<PathBuf>> {
    let mut inputs = args.input.clone();
    if let Some(list) = &args.input_list {
        inputs.extend(io::read_path_list(list)?);
    }
    if inputs.is_empty() {
        bail!("no input files: the input list is empty");
    }
    Ok(inputs)
}

/// Maps each input to `out_dir/<file name>.ksk`. Fails, before anything is written, if two
/// inputs would produce the same sketch file.
fn plan_jobs(inputs: Vec<PathBuf>, out_dir: &Path) -> Result<Vec<Job>> {
    let mut jobs = Vec::with_capacity(inputs.len());
    let mut by_output: HashMap<PathBuf, Vec<usize>> = HashMap::new();
    for input in inputs {
        let Some(file_name) = input.file_name() else {
            bail!("input {} is not a file name", input.display());
        };
        let mut name = file_name.to_owned();
        name.push(".");
        name.push(format::EXTENSION);
        let output = out_dir.join(name);
        by_output
            .entry(output.clone())
            .or_default()
            .push(jobs.len());
        jobs.push(Job { input, output });
    }
    let mut collisions: Vec<(PathBuf, Vec<usize>)> = by_output
        .into_iter()
        .filter(|(_, ix)| ix.len() > 1)
        .collect();
    if !collisions.is_empty() {
        collisions.sort();
        let shown: Vec<String> = collisions
            .iter()
            .take(5)
            .map(|(out, ix)| {
                let ins: Vec<String> = ix
                    .iter()
                    .map(|&i| jobs[i].input.display().to_string())
                    .collect();
                format!("  {} <- {}", out.display(), ins.join(", "))
            })
            .collect();
        let more = if collisions.len() > 5 { "\n  ..." } else { "" };
        bail!(
            "{} sketch file(s) would be written by more than one input. Rename the inputs or \
             sketch them into separate directories:\n{}{more}",
            collisions.len(),
            shown.join("\n"),
        );
    }
    Ok(jobs)
}

fn sketch_one(job: &Job, params: SketchParams) -> Result<Done> {
    let sketch = io::sketch_fastx(&job.input, params)?;
    if sketch.is_empty() {
        log::warn!(
            "{}: no k-mers sampled (empty file, or no contig of at least {} bp)",
            job.input.display(),
            params.min_contig_len,
        );
    }
    format::write(&sketch, &job.output)?;
    log::debug!(
        "{} -> {}: {} contigs ({} skipped), {} bp, {} sampled k-mers",
        job.input.display(),
        job.output.display(),
        sketch.num_contigs(),
        sketch.skipped_contigs(),
        sketch.total_length(),
        sketch.len(),
    );
    Ok(Done {
        contigs: sketch.num_contigs(),
        bases: sketch.total_length(),
        kmers: sketch.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn paths(v: &[&str]) -> Vec<PathBuf> {
        v.iter().map(PathBuf::from).collect()
    }

    #[test]
    fn outputs_keep_the_input_file_name() {
        let jobs = plan_jobs(paths(&["a.fa", "x/y/b.fna.gz"]), Path::new("out")).unwrap();
        assert_eq!(jobs[0].output, Path::new("out/a.fa.ksk"));
        assert_eq!(jobs[1].output, Path::new("out/b.fna.gz.ksk"));
    }

    #[test]
    fn same_file_name_in_two_directories_is_rejected() {
        let err = plan_jobs(paths(&["d1/a.fa", "d2/a.fa", "b.fa"]), Path::new("out"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("d1/a.fa") && err.contains("d2/a.fa"), "{err}");
    }

    #[test]
    fn duplicate_input_is_rejected() {
        assert!(plan_jobs(paths(&["a.fa", "a.fa"]), Path::new("out")).is_err());
    }
}
