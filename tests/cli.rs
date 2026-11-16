//! End-to-end checks of the `kadence` binary: exit codes and help text.

use std::process::{Command, Output};

fn kadence(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_kadence"))
        .args(args)
        .output()
        .expect("failed to run kadence")
}

fn stdout(o: &Output) -> String {
    String::from_utf8_lossy(&o.stdout).into_owned()
}

fn stderr(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).into_owned()
}

#[test]
fn help_lists_subcommands() {
    let o = kadence(&["--help"]);
    assert!(o.status.success());
    let out = stdout(&o);
    assert!(out.contains("sketch"), "{out}");
    assert!(out.contains("dist"), "{out}");
}

#[test]
fn version_matches_cargo() {
    let o = kadence(&["--version"]);
    assert!(o.status.success());
    assert_eq!(
        stdout(&o).trim(),
        format!("kadence {}", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn dist_help_documents_method_options_and_units() {
    let o = kadence(&["dist", "--help"]);
    assert!(o.status.success());
    let out = stdout(&o);
    for flag in [
        "--patience",
        "--min-anchors",
        "--break-alpha",
        "--min-ani",
        "--mode",
        "--kmer-size",
        "--scale",
        "--query-list",
        "--ref-list",
    ] {
        assert!(out.contains(flag), "missing {flag}:\n{out}");
    }
    // Help is wrapped to the terminal width, so compare with whitespace collapsed.
    let flat = out.split_whitespace().collect::<Vec<_>>().join(" ");
    assert!(
        flat.contains("non-matching sketched k-mers"),
        "patience units missing:\n{out}"
    );
}

#[test]
fn missing_arguments_exit_with_usage_error() {
    for args in [
        &[][..],
        &["sketch", "-o", "out"][..],
        &["dist", "-q", "a.fa"][..],
        &["dist", "-q", "a.fa", "-r", "b.fa", "--patience", "0"][..],
        &["dist", "-q", "a.fa", "-r", "b.fa", "--min-ani", "1.5"][..],
    ] {
        let o = kadence(args);
        assert_eq!(o.status.code(), Some(2), "{args:?}: {}", stderr(&o));
    }
}

#[test]
fn valid_invocations_report_not_implemented() {
    for (args, milestone) in [
        (&["sketch", "-i", "a.fa", "-o", "out"][..], "M2"),
        (&["dist", "-q", "a.fa", "-r", "b.fa", "-t", "4"][..], "M3"),
    ] {
        let o = kadence(args);
        assert_eq!(o.status.code(), Some(1), "{args:?}");
        let err = stderr(&o);
        assert!(err.contains("not implemented yet"), "{err}");
        assert!(err.contains(milestone), "{err}");
    }
}
