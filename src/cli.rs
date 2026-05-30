use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_multiinter::multiinter;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-multiinter", disable_help_flag = true)]
pub struct Cli {
    /// Input BED files (two or more, sorted by chrom then start).
    #[arg(short = 'i', long = "inputs", required = true, num_args = 1..)]
    pub inputs: Vec<PathBuf>,

    /// Names for each input file (same count as -i). Default: 1-based indices.
    #[arg(short = 'n', long = "names", num_args = 1..)]
    pub names: Vec<String>,

    /// Print a header line.
    #[arg(long = "header")]
    pub header: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let n = self.inputs.len();
        let names: Vec<String> = if self.names.is_empty() {
            (1..=n).map(|i| i.to_string()).collect()
        } else if self.names.len() != n {
            return Err(RsomicsError::Io(std::io::Error::other(format!(
                "--names count ({}) must match --inputs count ({})",
                self.names.len(),
                n
            ))));
        } else {
            self.names.clone()
        };

        let readers: Vec<File> = self
            .inputs
            .iter()
            .map(|p| File::open(p).map_err(RsomicsError::Io))
            .collect::<Result<Vec<_>>>()?;

        let stdout = io::stdout();
        let out = stdout.lock();
        multiinter(readers, &names, self.header, out)
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Report per-position intersection depth across N BED files (bedtools multiinter equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-i <FILE>... [OPTIONS]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "inputs",
                aliases: &[],
                value: Some("<path>..."),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BED files (sorted, two or more)",
                why_default: None,
            },
            FlagSpec {
                short: Some('n'),
                long: "names",
                aliases: &[],
                value: Some("<name>..."),
                type_hint: Some("String"),
                required: false,
                default: Some("1,2,...,N"),
                description: "Labels for each input file (must match count of -i)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "header",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Print column header line",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Find depth of coverage across three ChIP-seq peak files",
        command: "rsomics-bed-multiinter -i rep1.bed rep2.bed rep3.bed -n rep1 rep2 rep3",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
