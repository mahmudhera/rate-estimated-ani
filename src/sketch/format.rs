//! On-disk sketch format (`.ksk`), version 1.
//!
//! All integers are little-endian. Strings are a `u32` byte length followed by UTF-8 bytes.
//!
//! ```text
//! magic            8 bytes  "KDSKETCH"
//! version          u32      1
//! k                u8, then 3 zero bytes
//! c                u32
//! min_contig_len   u32
//! num_contigs      u32
//! num_entries      u32
//! skipped_contigs  u32
//! skipped_bases    u64
//! creator          string   e.g. "kadence 0.0.1"
//! name             string   input path as given
//! contigs          num_contigs x (name: string, length: u32)
//! contig_offsets   (num_contigs + 1) x u32
//! hashes           num_entries x u64      position order
//! positions        num_entries x u32      position order
//! forward          ceil(num_entries / 64) x u64, bit i = entry i
//! by_hash          num_entries x u32
//! checksum         u64      word-wise FNV-style hash of every preceding byte (see `checksum`)
//! ```
//!
//! Readers reject other versions, bad checksums, trailing bytes, and any sketch that fails
//! [`Sketch::validate`].

use std::fs;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

use super::{BitVec, Sketch};
use crate::params::SketchParams;

/// File signature.
pub const MAGIC: [u8; 8] = *b"KDSKETCH";
/// Current format version.
pub const FORMAT_VERSION: u32 = 1;
/// File extension for sketch files.
pub const EXTENSION: &str = "ksk";

const CREATOR: &str = concat!("kadence ", env!("CARGO_PKG_VERSION"));

/// Corruption check over the serialized bytes: FNV-1a constants applied to 8-byte
/// little-endian words (the tail is zero-padded), with a rotation so high bits also reach low
/// bits. Every step is a bijection of the state, so any single changed word changes the
/// result. It is about 8x faster than byte-wise FNV-1a. It is not a cryptographic hash.
fn checksum(bytes: &[u8]) -> u64 {
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut h: u64 = 0xcbf2_9ce4_8422_2325 ^ bytes.len() as u64;
    let mut words = bytes.chunks_exact(8);
    for w in &mut words {
        h = (h ^ u64::from_le_bytes(w.try_into().expect("8 bytes")))
            .wrapping_mul(PRIME)
            .rotate_left(29);
    }
    let tail = words.remainder();
    if !tail.is_empty() {
        let mut last = [0u8; 8];
        last[..tail.len()].copy_from_slice(tail);
        h = (h ^ u64::from_le_bytes(last))
            .wrapping_mul(PRIME)
            .rotate_left(29);
    }
    h
}

/// Serializes a sketch.
pub fn to_bytes(sketch: &Sketch) -> Vec<u8> {
    let n = sketch.len();
    let mut out = Vec::with_capacity(128 + n * 16 + sketch.name.len());
    let put_u32 = |out: &mut Vec<u8>, x: u32| out.extend_from_slice(&x.to_le_bytes());
    let put_str = |out: &mut Vec<u8>, s: &str| {
        out.extend_from_slice(&(s.len() as u32).to_le_bytes());
        out.extend_from_slice(s.as_bytes());
    };

    out.extend_from_slice(&MAGIC);
    put_u32(&mut out, FORMAT_VERSION);
    out.extend_from_slice(&[sketch.params.k, 0, 0, 0]);
    put_u32(&mut out, sketch.params.c);
    put_u32(&mut out, sketch.params.min_contig_len);
    put_u32(&mut out, sketch.num_contigs() as u32);
    put_u32(&mut out, n as u32);
    put_u32(&mut out, sketch.skipped_contigs);
    out.extend_from_slice(&sketch.skipped_bases.to_le_bytes());
    put_str(&mut out, CREATOR);
    put_str(&mut out, &sketch.name);
    for (name, &len) in sketch.contig_names.iter().zip(&sketch.contig_lengths) {
        put_str(&mut out, name);
        put_u32(&mut out, len);
    }
    for &o in &sketch.contig_offsets {
        put_u32(&mut out, o);
    }
    for &h in &sketch.hashes {
        out.extend_from_slice(&h.to_le_bytes());
    }
    for &p in &sketch.positions {
        put_u32(&mut out, p);
    }
    for &w in sketch.forward.words() {
        out.extend_from_slice(&w.to_le_bytes());
    }
    for &i in &sketch.by_hash {
        put_u32(&mut out, i);
    }
    let sum = checksum(&out);
    out.extend_from_slice(&sum.to_le_bytes());
    out
}

/// Bounds-checked little-endian reader.
struct Reader<'a> {
    buf: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.at.checked_add(n).filter(|&e| e <= self.buf.len());
        let Some(end) = end else {
            bail!("sketch file is truncated");
        };
        let s = &self.buf[self.at..end];
        self.at = end;
        Ok(s)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("4 bytes"),
        ))
    }

    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(
            self.take(8)?.try_into().expect("8 bytes"),
        ))
    }

    fn string(&mut self) -> Result<String> {
        let len = self.u32()? as usize;
        let bytes = self.take(len)?;
        String::from_utf8(bytes.to_vec()).context("string is not valid UTF-8")
    }

    /// Reads `n` u32s, checking the length first so a corrupt count can't force a huge
    /// allocation.
    fn u32s(&mut self, n: usize) -> Result<Vec<u32>> {
        let bytes = self.take(n.checked_mul(4).context("count overflow")?)?;
        Ok(bytes
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes(c.try_into().expect("4 bytes")))
            .collect())
    }

    fn u64s(&mut self, n: usize) -> Result<Vec<u64>> {
        let bytes = self.take(n.checked_mul(8).context("count overflow")?)?;
        Ok(bytes
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().expect("8 bytes")))
            .collect())
    }
}

/// Parses and validates a serialized sketch.
pub fn from_bytes(bytes: &[u8]) -> Result<Sketch> {
    ensure!(
        bytes.len() >= MAGIC.len() && bytes[..MAGIC.len()] == MAGIC,
        "not a kadence sketch (bad magic)"
    );
    ensure!(
        bytes.len() >= MAGIC.len() + 4 + 8,
        "sketch file is truncated"
    );
    let (body, tail) = bytes.split_at(bytes.len() - 8);
    let mut r = Reader {
        buf: body,
        at: MAGIC.len(),
    };
    let version = r.u32()?;
    ensure!(
        version == FORMAT_VERSION,
        "unsupported sketch format version {version} (this build reads version {FORMAT_VERSION})"
    );
    let stored = u64::from_le_bytes(tail.try_into().expect("8 bytes"));
    ensure!(
        stored == checksum(body),
        "sketch file is corrupt (checksum mismatch)"
    );

    let k = r.u8()?;
    ensure!(
        r.take(3)? == [0, 0, 0],
        "reserved header bytes are not zero"
    );
    let c = r.u32()?;
    let min_contig_len = r.u32()?;
    let num_contigs = r.u32()? as usize;
    let n = r.u32()? as usize;
    let skipped_contigs = r.u32()?;
    let skipped_bases = r.u64()?;
    let _creator = r.string()?;
    let name = r.string()?;
    // Each contig needs at least 8 bytes, so a corrupt count fails before allocating much.
    ensure!(
        num_contigs <= body.len() / 8,
        "contig count is larger than the file"
    );
    let mut contig_names = Vec::with_capacity(num_contigs);
    let mut contig_lengths = Vec::with_capacity(num_contigs);
    for _ in 0..num_contigs {
        contig_names.push(r.string()?);
        contig_lengths.push(r.u32()?);
    }
    let contig_offsets = r.u32s(num_contigs + 1)?;
    let hashes = r.u64s(n)?;
    let positions = r.u32s(n)?;
    let forward = BitVec::from_words(r.u64s(n.div_ceil(64))?, n)?;
    let by_hash = r.u32s(n)?;
    ensure!(r.at == body.len(), "unexpected bytes after the sketch data");

    let sketch = Sketch {
        name,
        params: SketchParams {
            k,
            c,
            min_contig_len,
        },
        contig_names,
        contig_lengths,
        contig_offsets,
        hashes,
        positions,
        forward,
        by_hash,
        skipped_contigs,
        skipped_bases,
    };
    sketch.validate().context("sketch file is inconsistent")?;
    Ok(sketch)
}

/// Writes a sketch to `path`. It writes to a temporary file in the same directory first and
/// then renames it, so a crash never leaves a partial sketch behind.
pub fn write(sketch: &Sketch, path: &Path) -> Result<()> {
    let bytes = to_bytes(sketch);
    let mut tmp_name = path.as_os_str().to_owned();
    tmp_name.push(".tmp");
    let tmp = Path::new(&tmp_name);
    let result = (|| -> Result<()> {
        // No fsync: renaming a fully written file is enough to never expose a partial sketch
        // after a crash of this process, and fsync per file is slow on network filesystems.
        let mut f = fs::File::create(tmp)?;
        f.write_all(&bytes)?;
        drop(f);
        fs::rename(tmp, path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(tmp);
    }
    result.with_context(|| format!("failed to write sketch {}", path.display()))
}

/// Reads and validates the sketch at `path`.
pub fn read(path: &Path) -> Result<Sketch> {
    let bytes = fs::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    from_bytes(&bytes).with_context(|| format!("failed to load sketch {}", path.display()))
}

/// Whether `path` starts with the sketch signature. Used to tell sketches from FASTA input.
pub fn is_sketch_file(path: &Path) -> Result<bool> {
    use std::io::Read;
    let mut f =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut head = [0u8; MAGIC.len()];
    match f.read_exact(&mut head) {
        Ok(()) => Ok(head == MAGIC),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e).with_context(|| format!("failed to read {}", path.display())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::SketchBuilder;
    use crate::sketch::kmer::tests::random_seq;

    fn sample() -> Sketch {
        let mut b = SketchBuilder::new(
            "dir/genome.fa",
            SketchParams {
                k: 17,
                c: 20,
                min_contig_len: 100,
            },
        );
        b.add_contig("ctg1 some description", &random_seq(30_000, 1))
            .unwrap();
        b.add_contig("short", b"ACGTACGT").unwrap();
        b.add_contig("ctg2", &random_seq(12_345, 2)).unwrap();
        b.finish()
    }

    #[test]
    fn round_trip() {
        let s = sample();
        let back = from_bytes(&to_bytes(&s)).unwrap();
        assert_eq!(back, s);
    }

    #[test]
    fn round_trip_empty() {
        let s = SketchBuilder::new("empty.fa", SketchParams::default()).finish();
        assert_eq!(from_bytes(&to_bytes(&s)).unwrap(), s);
    }

    #[test]
    fn round_trip_through_a_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("g.ksk");
        let s = sample();
        write(&s, &path).unwrap();
        assert!(is_sketch_file(&path).unwrap());
        assert_eq!(read(&path).unwrap(), s);
        // No temporary file left behind.
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn rejects_corruption() {
        let bytes = to_bytes(&sample());
        // Flipped bit in the payload.
        let mut bad = bytes.clone();
        let mid = bad.len() / 2;
        bad[mid] ^= 1;
        let err = from_bytes(&bad).unwrap_err().to_string();
        assert!(err.contains("checksum"), "{err}");
        // Truncation.
        assert!(from_bytes(&bytes[..bytes.len() - 20]).is_err());
        // Wrong magic.
        let mut bad = bytes.clone();
        bad[0] = b'X';
        assert!(from_bytes(&bad).unwrap_err().to_string().contains("magic"));
        // Future version.
        let mut bad = bytes.clone();
        bad[8..12].copy_from_slice(&2u32.to_le_bytes());
        assert!(
            from_bytes(&bad)
                .unwrap_err()
                .to_string()
                .contains("version")
        );
        // Trailing garbage.
        let mut bad = bytes;
        bad.push(0);
        assert!(from_bytes(&bad).is_err());
    }

    #[test]
    fn fasta_is_not_a_sketch() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("x.fa");
        fs::write(&fa, ">a\nACGT\n").unwrap();
        assert!(!is_sketch_file(&fa).unwrap());
        let tiny = dir.path().join("tiny");
        fs::write(&tiny, "KD").unwrap();
        assert!(!is_sketch_file(&tiny).unwrap());
    }
}
