//! Multi-file interval intersection depth sweep.
//!
//! ## Algorithm
//!
//! For each position in the genome, reports how many of the N input BED files
//! have at least one interval covering that position. Uses a coordinate-sweep
//! (priority-queue over start/end events) per chromosome, emitting a row for
//! every contiguous run of constant depth.
//!
//! Output columns (tab-separated):
//! ```text
//! chrom  start  end  depth  list  file1  file2  ...  fileN
//! ```
//! - `depth`: number of files with coverage at this position.
//! - `list`: comma-separated names (1-based indices, or `-names` labels) of
//!   covering files; `"none"` when depth = 0.
//! - `file1..fileN`: 1 if that file has coverage, 0 otherwise.
//!
//! All input files must be sorted by chrom then start. Chromosomes are
//! processed in the order they first appear across files (matching bedtools).
//!
//! ## Reference
//!
//! `BEDTools` multiinter — Quinlan & Hall (2010). Bioinformatics 26(6): 841–842.
//! DOI: 10.1093/bioinformatics/btq033

use std::collections::BinaryHeap;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use rsomics_common::{Result, RsomicsError};

// ── Data types ─────────────────────────────────────────────────────────────

/// A parsed BED record (first three columns only).
#[derive(Debug, Clone)]
struct Bed3 {
    chrom: String,
    start: i64,
    end: i64,
}

/// One coordinate event: a start or end of a file's interval.
#[derive(Debug, Eq, PartialEq)]
struct Event {
    /// Coordinate position.
    coord: i64,
    /// Whether this is a start (true) or end (false) event.
    is_start: bool,
    /// Which file (0-based index) this event belongs to.
    file_idx: usize,
}

/// Events are ordered by coord ascending, then ends before starts at the same
/// coord (matching bedtools' tie-breaking: process END before START so an
/// interval [10,20) and [20,30) do not share a reported segment at coord 20).
impl Ord for Event {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // BinaryHeap is a max-heap; we want min-heap, so reverse.
        other
            .coord
            .cmp(&self.coord)
            .then(self.is_start.cmp(&other.is_start)) // false (end) < true (start)
    }
}
impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// ── Reader ─────────────────────────────────────────────────────────────────

struct FileReader<R: BufRead> {
    reader: R,
    /// The buffered "next" record (peeked but not yet consumed).
    peeked: Option<Bed3>,
}

impl<R: BufRead> FileReader<R> {
    fn new(reader: R) -> Result<Self> {
        let mut fr = FileReader {
            reader,
            peeked: None,
        };
        fr.peeked = fr.read_next()?;
        Ok(fr)
    }

    fn read_next(&mut self) -> Result<Option<Bed3>> {
        loop {
            let mut line = String::new();
            let n = self.reader.read_line(&mut line).map_err(RsomicsError::Io)?;
            if n == 0 {
                return Ok(None);
            }
            let t = line.trim_end_matches(['\n', '\r']);
            if t.is_empty()
                || t.starts_with('#')
                || t.starts_with("track")
                || t.starts_with("browser")
            {
                continue;
            }
            let mut cols = t.splitn(4, '\t');
            let chrom = match cols.next() {
                Some(c) => c.to_owned(),
                None => continue,
            };
            let start: i64 = match cols.next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let end: i64 = match cols.next().and_then(|s| s.parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            return Ok(Some(Bed3 { chrom, start, end }));
        }
    }

    /// Peek at the current record without consuming it.
    fn peek(&self) -> Option<&Bed3> {
        self.peeked.as_ref()
    }

    /// Consume and return the current record, loading the next.
    fn next(&mut self) -> Result<Option<Bed3>> {
        let cur = self.peeked.take();
        self.peeked = self.read_next()?;
        Ok(cur)
    }
}

// ── Core sweep ─────────────────────────────────────────────────────────────

/// Run the multi-interval intersection sweep.
///
/// Reads N BED files from `readers` (each `(impl Read)`), sweeps coordinates,
/// and writes the result table to `out`.
pub fn multiinter<R: Read, W: Write>(
    readers: Vec<R>,
    names: &[String],
    header: bool,
    out: W,
) -> Result<()> {
    let n = readers.len();
    let mut file_readers: Vec<FileReader<BufReader<R>>> = readers
        .into_iter()
        .map(|r| FileReader::new(BufReader::new(r)))
        .collect::<Result<Vec<_>>>()?;

    let mut writer = BufWriter::new(out);

    if header {
        write!(writer, "chrom\tstart\tend\tnum\tlist").map_err(RsomicsError::Io)?;
        for name in names {
            write!(writer, "\t{name}").map_err(RsomicsError::Io)?;
        }
        writeln!(writer).map_err(RsomicsError::Io)?;
    }

    // Determine chromosome processing order: first appearance across all files.
    // We collect the first chrom of each file (their "current" chrom) and build
    // an ordered list as we process.
    let mut processed_chroms: Vec<String> = Vec::new();

    loop {
        // Find the lexicographically smallest chrom among all files' current records.
        let next_chrom = file_readers
            .iter()
            .filter_map(|fr| fr.peek().map(|r| r.chrom.as_str()))
            .min()
            .map(ToOwned::to_owned);

        let Some(chrom) = next_chrom else {
            break;
        };

        if !processed_chroms.contains(&chrom) {
            processed_chroms.push(chrom.clone());
        }

        // Build the event queue for this chromosome from all files.
        let mut queue: BinaryHeap<Event> = BinaryHeap::new();

        for (file_idx, fr) in file_readers.iter_mut().enumerate() {
            // Consume all intervals on this chromosome, pushing events.
            while fr.peek().is_some_and(|r| r.chrom == chrom) {
                let rec = fr.next()?.unwrap();
                queue.push(Event {
                    coord: rec.start,
                    is_start: true,
                    file_idx,
                });
                queue.push(Event {
                    coord: rec.end,
                    is_start: false,
                    file_idx,
                });
            }
        }

        // Sweep: track which files are "active" (have depth > 0) at each coord.
        // Per-file depth tracks nesting so we know exactly when presence changes.
        // A new output row is emitted only when the PRESENCE SET changes — matching
        // bedtools' behaviour of outputting 0/1 per file, not raw depth.
        let mut depth = vec![0u32; n];
        let mut active_count = 0usize;
        let mut segment_start: i64 = 0;
        let mut in_segment = false;

        while let Some(event) = queue.pop() {
            let coord = event.coord;

            // Does this event change which files are PRESENT (depth 0→1 or 1→0)?
            let presence_changes = if event.is_start {
                depth[event.file_idx] == 0
            } else {
                depth[event.file_idx] == 1
            };

            // When presence set changes: close the current open segment (if any).
            if presence_changes && in_segment && coord > segment_start {
                emit_row(
                    &mut writer,
                    &chrom,
                    segment_start,
                    coord,
                    active_count,
                    &depth,
                    names,
                )?;
            }

            // Update depth.
            if event.is_start {
                depth[event.file_idx] += 1;
                if depth[event.file_idx] == 1 {
                    active_count += 1;
                }
            } else {
                depth[event.file_idx] -= 1;
                if depth[event.file_idx] == 0 {
                    active_count = active_count.saturating_sub(1);
                }
            }

            // After updating, start a new segment at coord if there are active files.
            if presence_changes {
                if active_count > 0 {
                    segment_start = coord;
                    in_segment = true;
                } else {
                    in_segment = false;
                }
            }
        }
    }

    writer.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

fn emit_row<W: Write>(
    w: &mut W,
    chrom: &str,
    start: i64,
    end: i64,
    active_count: usize,
    depth: &[u32],
    names: &[String],
) -> Result<()> {
    // Build the comma-separated list of covering file names.
    let mut list = String::new();
    let mut first = true;
    for (i, &d) in depth.iter().enumerate() {
        if d > 0 {
            if !first {
                list.push(',');
            }
            list.push_str(&names[i]);
            first = false;
        }
    }
    if list.is_empty() {
        list.push_str("none");
    }

    write!(w, "{chrom}\t{start}\t{end}\t{active_count}\t{list}").map_err(RsomicsError::Io)?;
    for &d in depth {
        write!(w, "\t{}", u32::from(d > 0)).map_err(RsomicsError::Io)?;
    }
    writeln!(w).map_err(RsomicsError::Io)
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn run(inputs: &[&str], names: &[&str]) -> String {
        let readers: Vec<Cursor<&str>> = inputs.iter().map(|s| Cursor::new(*s)).collect();
        let name_strings: Vec<String> = names.iter().copied().map(str::to_string).collect();
        let mut out = Vec::new();
        multiinter(readers, &name_strings, false, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn no_overlap_two_files() {
        let out = run(&["chr1\t100\t200\n", "chr1\t300\t400\n"], &["f1", "f2"]);
        let lines: Vec<&str> = out.lines().collect();
        // Two non-overlapping intervals → 2 rows, each with depth 1
        assert_eq!(lines.len(), 2, "lines: {out:?}");
        assert!(lines[0].contains("\t1\tf1\t1\t0"), "row0: {}", lines[0]);
        assert!(lines[1].contains("\t1\tf2\t0\t1"), "row1: {}", lines[1]);
    }

    #[test]
    fn full_overlap_two_files() {
        let out = run(&["chr1\t100\t300\n", "chr1\t100\t300\n"], &["f1", "f2"]);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 1, "lines: {out:?}");
        assert!(lines[0].contains("\t2\tf1,f2\t1\t1"), "row: {}", lines[0]);
    }

    #[test]
    fn partial_overlap_three_files() {
        // A:[100,200) B:[150,250) C:[50,120)
        // Events: 50(S,C),100(S,A),120(E,C),150(S,B),200(E,A),250(E,B)
        // Segments: [50,100)→C only, [100,120)→A+C, [120,150)→A only,
        //           [150,200)→A+B, [200,250)→B only
        let out = run(
            &["chr1\t100\t200\n", "chr1\t150\t250\n", "chr1\t50\t120\n"],
            &["f1", "f2", "f3"],
        );
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 5, "lines: {out:?}");
        // [50,100) depth=1, list=f3
        assert!(lines[0].starts_with("chr1\t50\t100\t1\tf3"), "{}", lines[0]);
        // [100,120) depth=2, list=f1,f3
        assert!(
            lines[1].starts_with("chr1\t100\t120\t2\tf1,f3"),
            "{}",
            lines[1]
        );
        // [120,150) depth=1, list=f1
        assert!(
            lines[2].starts_with("chr1\t120\t150\t1\tf1"),
            "{}",
            lines[2]
        );
        // [150,200) depth=2, list=f1,f2
        assert!(
            lines[3].starts_with("chr1\t150\t200\t2\tf1,f2"),
            "{}",
            lines[3]
        );
        // [200,250) depth=1, list=f2
        assert!(
            lines[4].starts_with("chr1\t200\t250\t1\tf2"),
            "{}",
            lines[4]
        );
    }

    #[test]
    fn header_printed_when_requested() {
        let readers: Vec<Cursor<&str>> = vec![Cursor::new("chr1\t10\t20\n")];
        let names = vec!["s1".to_string()];
        let mut out = Vec::new();
        multiinter(readers, &names, true, &mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.starts_with("chrom\tstart\tend\tnum\tlist\ts1"), "{s}");
    }

    #[test]
    fn multi_chrom_ordering() {
        let out = run(
            &["chr1\t10\t20\nchr2\t10\t20\n", "chr2\t15\t25\n"],
            &["f1", "f2"],
        );
        let lines: Vec<&str> = out.lines().collect();
        // chr1: f1 only → 1 row
        // chr2: [10,15) f1, [15,20) f1+f2, [20,25) f2 → 3 rows
        assert_eq!(lines.len(), 4, "lines:\n{out}");
        assert!(lines[0].starts_with("chr1"), "{}", lines[0]);
        assert!(lines[1].starts_with("chr2\t10\t15"), "{}", lines[1]);
        assert!(lines[2].starts_with("chr2\t15\t20"), "{}", lines[2]);
        assert!(lines[3].starts_with("chr2\t20\t25"), "{}", lines[3]);
    }

    #[test]
    fn empty_file_does_not_crash() {
        let out = run(&["", "chr1\t100\t200\n"], &["f1", "f2"]);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 1);
        assert!(lines[0].contains("\t1\tf2\t0\t1"), "{}", lines[0]);
    }
}
