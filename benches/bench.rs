use std::io::Cursor;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_bed_multiinter::multiinter;

fn make_bed(n: usize, seed: u64) -> String {
    use std::fmt::Write as FmtWrite;
    let mut rng = seed;
    let mut out = String::with_capacity(n * 30);
    let chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"];
    let chrom_size = 200_000_000i64;

    let mut records: Vec<(usize, i64, i64)> = (0..n)
        .map(|_| {
            rng = rng
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let ci = ((rng >> 33) as usize) % chroms.len();
            rng = rng
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let start = ((rng >> 33) as i64) % (chrom_size - 2000);
            let len = 100 + ((rng & 0xFFF) as i64);
            (ci, start, start + len)
        })
        .collect();
    records.sort_by_key(|r| (r.0, r.1));
    for (ci, s, e) in records {
        let _ = writeln!(out, "{}\t{}\t{}", chroms[ci], s, e);
    }
    out
}

fn bench_multiinter(c: &mut Criterion) {
    let bed1 = make_bed(200_000, 42);
    let bed2 = make_bed(200_000, 137);
    let bed3 = make_bed(200_000, 9999);
    let names: Vec<String> = vec!["f1".to_string(), "f2".to_string(), "f3".to_string()];

    c.bench_function("multiinter_3x200k", |b| {
        b.iter(|| {
            let readers = vec![
                Cursor::new(bed1.as_str()),
                Cursor::new(bed2.as_str()),
                Cursor::new(bed3.as_str()),
            ];
            let mut out = Vec::new();
            multiinter(readers, &names, false, &mut out).unwrap();
        });
    });
}

criterion_group!(benches, bench_multiinter);
criterion_main!(benches);
