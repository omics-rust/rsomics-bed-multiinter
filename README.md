# rsomics-bed-multiinter

Report per-position intersection depth across N sorted BED files — `bedtools multiinter` equivalent.

## Usage

```
rsomics-bed-multiinter -i rep1.bed rep2.bed rep3.bed [-n rep1 rep2 rep3] [--header]
```

## Output

Tab-separated: `chrom start end depth list file1 file2 ... fileN`

- `depth`: number of files covering this position
- `list`: comma-separated names of covering files (or `none`)
- `file1..fileN`: 1/0 per file

## Origin

Independent Rust reimplementation of `bedtools multiinter` based on the
published method and format specification.

License: MIT OR Apache-2.0  
Upstream credit: [BEDTools](https://github.com/arq5x/bedtools2) (MIT / GPL-2.0)
