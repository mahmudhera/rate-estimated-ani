//! End-to-end tests of `kadence sketch`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use kadence::sketch::format;

fn kadence(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_kadence"))
        .args(args)
        .output()
        .expect("failed to run kadence")
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

/// Deterministic pseudo-random ACGT sequence.
fn random_seq(len: usize, mut state: u64) -> String {
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ['A', 'C', 'G', 'T'][(state >> 32) as usize & 3]
        })
        .collect()
}

fn write_fasta(path: &Path, contigs: &[(&str, usize, u64)]) {
    let mut text = String::new();
    for (name, len, seed) in contigs {
        text.push_str(&format!(">{name}\n{}\n", random_seq(*len, *seed)));
    }
    fs::write(path, text).unwrap();
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

#[test]
fn sketches_files_and_list_inputs() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.fa");
    let b = dir.path().join("b.fna");
    let c = dir.path().join("c.fa");
    write_fasta(&a, &[("a1", 20_000, 1), ("a2", 300, 2)]);
    write_fasta(&b, &[("b1", 15_000, 3)]);
    write_fasta(&c, &[("c1", 9_000, 4)]);
    let list = dir.path().join("list.txt");
    fs::write(&list, format!("{}\n\n{}\n", p(&b), p(&c))).unwrap();
    let out = dir.path().join("sketches");

    let o = kadence(&[
        "sketch",
        "-i",
        p(&a),
        "--input-list",
        p(&list),
        "-o",
        p(&out),
        "-k",
        "19",
        "-c",
        "50",
        "-t",
        "2",
    ]);
    assert!(o.status.success(), "{}", stderr(&o));

    let sa = format::read(&out.join("a.fa.ksk")).unwrap();
    assert_eq!(sa.name(), p(&a));
    assert_eq!((sa.params().k, sa.params().c), (19, 50));
    // The 300 bp contig is under the default 500 bp minimum.
    assert_eq!(sa.contig_names(), ["a1"]);
    assert_eq!(sa.total_length(), 20_000);
    assert_eq!(sa.skipped_contigs(), 1);
    assert!(!sa.is_empty());
    for name in ["b.fna.ksk", "c.fa.ksk"] {
        assert!(format::read(&out.join(name)).is_ok(), "{name}");
    }
}

#[test]
fn min_contig_len_flag_keeps_short_contigs() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.fa");
    write_fasta(&a, &[("long", 5_000, 1), ("short", 300, 2)]);
    let out = dir.path().join("out");
    let o = kadence(&[
        "sketch",
        "-i",
        p(&a),
        "-o",
        p(&out),
        "--min-contig-len",
        "0",
    ]);
    assert!(o.status.success(), "{}", stderr(&o));
    let s = format::read(&out.join("a.fa.ksk")).unwrap();
    assert_eq!(s.num_contigs(), 2);
    assert_eq!(s.params().min_contig_len, 0);
}

#[test]
fn output_is_identical_across_thread_counts() {
    let dir = tempfile::tempdir().unwrap();
    let inputs: Vec<PathBuf> = (0..6)
        .map(|i| {
            let f = dir.path().join(format!("g{i}.fa"));
            write_fasta(&f, &[("x", 10_000 + i * 1_000, i as u64 + 1)]);
            f
        })
        .collect();
    let mut outs = Vec::new();
    for threads in ["1", "4"] {
        let out = dir.path().join(format!("out{threads}"));
        let mut args = vec!["sketch", "-o", p(&out), "-t", threads, "-i"];
        args.extend(inputs.iter().map(|f| p(f)));
        let o = kadence(&args);
        assert!(o.status.success(), "{}", stderr(&o));
        outs.push(out);
    }
    for i in 0..6 {
        let name = format!("g{i}.fa.ksk");
        assert_eq!(
            fs::read(outs[0].join(&name)).unwrap(),
            fs::read(outs[1].join(&name)).unwrap(),
            "{name}"
        );
    }
}

#[test]
fn name_collision_fails_before_writing() {
    let dir = tempfile::tempdir().unwrap();
    for sub in ["d1", "d2"] {
        fs::create_dir(dir.path().join(sub)).unwrap();
        write_fasta(&dir.path().join(sub).join("a.fa"), &[("x", 2_000, 1)]);
    }
    let out = dir.path().join("out");
    let o = kadence(&[
        "sketch",
        "-i",
        p(&dir.path().join("d1/a.fa")),
        p(&dir.path().join("d2/a.fa")),
        "-o",
        p(&out),
    ]);
    assert_eq!(o.status.code(), Some(1));
    assert!(stderr(&o).contains("more than one input"), "{}", stderr(&o));
    assert!(!out.exists(), "nothing should be written");
}

#[test]
fn bad_inputs_fail_but_good_ones_are_written() {
    let dir = tempfile::tempdir().unwrap();
    let good = dir.path().join("good.fa");
    write_fasta(&good, &[("x", 3_000, 1)]);
    let bad = dir.path().join("bad.fa");
    fs::write(&bad, "not a fasta file\n").unwrap();
    let missing = dir.path().join("missing.fa");
    let out = dir.path().join("out");
    let o = kadence(&[
        "sketch",
        "-i",
        p(&good),
        p(&bad),
        p(&missing),
        "-o",
        p(&out),
    ]);
    assert_eq!(o.status.code(), Some(1));
    let err = stderr(&o);
    assert!(err.contains("2 of 3 input file(s) failed"), "{err}");
    assert!(
        err.contains("bad.fa") && err.contains("missing.fa"),
        "{err}"
    );
    assert!(out.join("good.fa.ksk").exists());
    assert!(!out.join("bad.fa.ksk").exists());
}

#[test]
fn existing_sketches_are_overwritten() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.fa");
    write_fasta(&a, &[("x", 4_000, 1)]);
    let out = dir.path().join("out");
    fs::create_dir(&out).unwrap();
    fs::write(out.join("a.fa.ksk"), "stale").unwrap();
    let o = kadence(&["sketch", "-i", p(&a), "-o", p(&out)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(format::read(&out.join("a.fa.ksk")).is_ok());
}

#[test]
fn empty_input_warns_and_writes_an_empty_sketch() {
    let dir = tempfile::tempdir().unwrap();
    let e = dir.path().join("empty.fa");
    fs::write(&e, "").unwrap();
    let out = dir.path().join("out");
    let o = kadence(&["sketch", "-i", p(&e), "-o", p(&out)]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(stderr(&o).contains("no k-mers sampled"), "{}", stderr(&o));
    assert!(format::read(&out.join("empty.fa.ksk")).unwrap().is_empty());
}
