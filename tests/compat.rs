/// Compatibility tests for rsomics-bed-multiinter vs bedtools multiinter.
///
/// Each test writes small BED fixtures to a temp dir, runs both binaries, and
/// compares tab-separated output field by field.
use std::process::Command;

fn bedtools_ok() -> bool {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn run_ours(args: &[&str]) -> String {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-multiinter");
    let out = Command::new(bin)
        .args(args)
        .output()
        .expect("rsomics-bed-multiinter failed");
    assert!(
        out.status.success(),
        "exit {:?}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn run_bedtools(args: &[&str]) -> String {
    let out = Command::new("bedtools")
        .args(args)
        .output()
        .expect("bedtools failed");
    assert!(
        out.status.success(),
        "bedtools exit {:?}\nstderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn write_tmp(name: &str, content: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir();
    let p = dir.join(name);
    std::fs::write(&p, content).unwrap();
    p
}

/// Parse multiinter output into (chrom, start, end, depth) per row, skipping header.
fn parse_table(s: &str) -> Vec<(String, i64, i64, usize)> {
    s.lines()
        .filter(|l| !l.starts_with("chrom"))
        .filter(|l| !l.is_empty())
        .map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            assert!(cols.len() >= 4, "bad row: {l}");
            (
                cols[0].to_owned(),
                cols[1].parse().unwrap(),
                cols[2].parse().unwrap(),
                cols[3].parse().unwrap(),
            )
        })
        .collect()
}

#[test]
fn no_overlap_two_files() {
    let a = write_tmp("multiinter_compat_a1.bed", "chr1\t100\t200\n");
    let b = write_tmp("multiinter_compat_b1.bed", "chr1\t300\t400\n");

    let ours = run_ours(&["-i", a.to_str().unwrap(), b.to_str().unwrap()]);
    let rows = parse_table(&ours);
    assert_eq!(rows.len(), 2, "expected 2 non-overlapping rows: {ours:?}");
    assert_eq!(rows[0].3, 1, "depth of first row");
    assert_eq!(rows[1].3, 1, "depth of second row");
}

#[test]
fn full_overlap_two_files() {
    let a = write_tmp("multiinter_compat_a2.bed", "chr1\t100\t300\n");
    let b = write_tmp("multiinter_compat_b2.bed", "chr1\t100\t300\n");

    let ours = run_ours(&["-i", a.to_str().unwrap(), b.to_str().unwrap()]);
    let rows = parse_table(&ours);
    assert_eq!(rows.len(), 1, "full overlap → single row: {ours:?}");
    assert_eq!(rows[0].3, 2, "depth = 2 when both files cover same region");
}

#[test]
fn column_count_matches_bedtools() {
    if !bedtools_ok() {
        eprintln!("bedtools not found, skipping");
        return;
    }
    let a = write_tmp("multiinter_compat_a3.bed", "chr1\t10\t50\nchr1\t60\t90\n");
    let b = write_tmp("multiinter_compat_b3.bed", "chr1\t30\t70\n");

    let ours_out = run_ours(&["-i", a.to_str().unwrap(), b.to_str().unwrap()]);
    let bt_out = run_bedtools(&["multiinter", "-i", a.to_str().unwrap(), b.to_str().unwrap()]);

    let ours_rows = parse_table(&ours_out);
    let bt_rows = parse_table(&bt_out);
    assert_eq!(
        ours_rows.len(),
        bt_rows.len(),
        "row counts differ.\nours:\n{ours_out}\nbedtools:\n{bt_out}"
    );
}

#[test]
fn depth_matches_bedtools_per_row() {
    if !bedtools_ok() {
        eprintln!("bedtools not found, skipping");
        return;
    }
    let a = write_tmp(
        "multiinter_compat_a4.bed",
        "chr1\t100\t200\nchr2\t50\t150\n",
    );
    let b = write_tmp(
        "multiinter_compat_b4.bed",
        "chr1\t150\t250\nchr2\t50\t100\n",
    );
    let c = write_tmp("multiinter_compat_c4.bed", "chr1\t10\t300\n");

    let ours_out = run_ours(&[
        "-i",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        c.to_str().unwrap(),
    ]);
    let bt_out = run_bedtools(&[
        "multiinter",
        "-i",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        c.to_str().unwrap(),
    ]);

    let ours_rows = parse_table(&ours_out);
    let bt_rows = parse_table(&bt_out);
    assert_eq!(
        ours_rows.len(),
        bt_rows.len(),
        "row count mismatch.\nours:\n{ours_out}\nbedtools:\n{bt_out}"
    );

    for (i, (o, b)) in ours_rows.iter().zip(bt_rows.iter()).enumerate() {
        assert_eq!(o.0, b.0, "row {i} chrom mismatch");
        assert_eq!(o.1, b.1, "row {i} start mismatch");
        assert_eq!(o.2, b.2, "row {i} end mismatch");
        assert_eq!(o.3, b.3, "row {i} depth mismatch");
    }
}

#[test]
fn header_flag_adds_header_line() {
    let a = write_tmp("multiinter_compat_a5.bed", "chr1\t10\t20\n");
    let b = write_tmp("multiinter_compat_b5.bed", "chr1\t15\t25\n");

    let with_header = run_ours(&["-i", a.to_str().unwrap(), b.to_str().unwrap(), "--header"]);
    let without_header = run_ours(&["-i", a.to_str().unwrap(), b.to_str().unwrap()]);

    assert!(
        with_header.starts_with("chrom\t"),
        "header not present: {with_header:?}"
    );
    assert!(
        !without_header.starts_with("chrom\t"),
        "unexpected header: {without_header:?}"
    );

    let header_lines: Vec<&str> = with_header.lines().collect();
    let no_header_lines: Vec<&str> = without_header.lines().collect();
    assert_eq!(
        header_lines.len(),
        no_header_lines.len() + 1,
        "with_header should have exactly one extra line"
    );
}
