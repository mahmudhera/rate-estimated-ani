//! Reading inputs: path lists and FASTA/FASTQ files.
//!
//! Sequence parsing follows skani v0.4.0 (`src/file_io.rs`; MIT OR Apache-2.0, Copyright (c)
//! 2022 Jim Shaw): needletail reads plain or compressed FASTA/FASTQ, and each record becomes
//! one contig of the genome.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use needletail::errors::ParseErrorKind;

use crate::params::SketchParams;
use crate::sketch::{Sketch, SketchBuilder};

/// Reads a file that lists one path per line. Surrounding whitespace is trimmed and blank
/// lines are skipped. Relative paths are kept as written (relative to the working directory).
pub fn read_path_list(path: &Path) -> Result<Vec<PathBuf>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read list file {}", path.display()))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(PathBuf::from)
        .collect())
}

/// Sketches a FASTA/FASTQ file (plain, gzip, bzip2, xz or zstd). Every record is a contig. An
/// empty file gives an empty sketch.
pub fn sketch_fastx(path: &Path, params: SketchParams) -> Result<Sketch> {
    let mut builder = SketchBuilder::new(path.to_string_lossy(), params);
    let mut reader = match needletail::parse_fastx_file(path) {
        Ok(r) => r,
        Err(e) if e.kind == ParseErrorKind::EmptyFile => return Ok(builder.finish()),
        Err(e) => {
            return Err(e).with_context(|| format!("failed to open {}", path.display()));
        }
    };
    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        let name = String::from_utf8_lossy(record.id());
        builder
            .add_contig(&name, &record.seq())
            .with_context(|| format!("failed to sketch {}", path.display()))?;
    }
    Ok(builder.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::kmer::tests::random_seq;
    use std::io::Write;

    fn params(min_contig_len: u32) -> SketchParams {
        SketchParams {
            k: 17,
            c: 10,
            min_contig_len,
        }
    }

    fn fasta(records: &[(&str, &[u8])]) -> Vec<u8> {
        let mut out = Vec::new();
        for (name, seq) in records {
            out.extend_from_slice(format!(">{name}\n").as_bytes());
            // Wrap at 60 columns, as most FASTA files do.
            for line in seq.chunks(60) {
                out.extend_from_slice(line);
                out.push(b'\n');
            }
        }
        out
    }

    #[test]
    fn path_list_skips_blank_lines() {
        let dir = tempfile::tempdir().unwrap();
        let list = dir.path().join("list.txt");
        fs::write(&list, "a.fa\n\n  dir/b.fa  \n\t\nc.fna.gz").unwrap();
        assert_eq!(
            read_path_list(&list).unwrap(),
            ["a.fa", "dir/b.fa", "c.fna.gz"].map(PathBuf::from)
        );
    }

    #[test]
    fn multi_record_fasta_matches_direct_sketch() {
        let a = random_seq(5_000, 1);
        let b = random_seq(3_000, 2);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("g.fa");
        fs::write(
            &path,
            fasta(&[("a desc", &a), ("tiny", b"ACGT"), ("b", &b)]),
        )
        .unwrap();

        let s = sketch_fastx(&path, params(100)).unwrap();
        let mut direct = SketchBuilder::new(path.to_string_lossy(), params(100));
        direct.add_contig("a desc", &a).unwrap();
        direct.add_contig("tiny", b"ACGT").unwrap();
        direct.add_contig("b", &b).unwrap();
        assert_eq!(s, direct.finish());
        assert_eq!(s.contig_names(), ["a desc", "b"]);
        assert_eq!(s.skipped_contigs(), 1);
    }

    #[test]
    fn gzip_input_matches_plain() {
        let seq = random_seq(8_000, 3);
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("g.fa");
        let gz = dir.path().join("g.fa.gz");
        let data = fasta(&[("x", &seq)]);
        fs::write(&plain, &data).unwrap();
        let mut enc = flate2::write::GzEncoder::new(
            fs::File::create(&gz).unwrap(),
            flate2::Compression::default(),
        );
        enc.write_all(&data).unwrap();
        enc.finish().unwrap();

        let a = sketch_fastx(&plain, params(0)).unwrap();
        let b = sketch_fastx(&gz, params(0)).unwrap();
        assert_eq!(a.hashes(), b.hashes());
        assert_eq!(a.positions(), b.positions());
    }

    #[test]
    fn fastq_records_are_contigs() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.fq");
        let seq = random_seq(200, 4);
        let qual = vec![b'I'; seq.len()];
        let mut data = b"@r1\n".to_vec();
        data.extend_from_slice(&seq);
        data.extend_from_slice(b"\n+\n");
        data.extend_from_slice(&qual);
        data.push(b'\n');
        fs::write(&path, data).unwrap();
        let s = sketch_fastx(&path, params(0)).unwrap();
        assert_eq!(s.num_contigs(), 1);
        assert_eq!(s.total_length(), 200);
    }

    #[test]
    fn empty_file_gives_empty_sketch() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.fa");
        fs::write(&path, "").unwrap();
        let s = sketch_fastx(&path, params(0)).unwrap();
        assert!(s.is_empty());
        assert_eq!(s.num_contigs(), 0);
    }

    #[test]
    fn missing_and_malformed_files_are_errors() {
        let dir = tempfile::tempdir().unwrap();
        assert!(sketch_fastx(&dir.path().join("nope.fa"), params(0)).is_err());
        let bad = dir.path().join("bad.fa");
        fs::write(&bad, "this is not fasta\n").unwrap();
        assert!(sketch_fastx(&bad, params(0)).is_err());
    }
}
