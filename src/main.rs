use anyhow::Result;
use clap::{Args, Parser, Subcommand};

mod bedcov;
mod calibrate;
mod region;

#[derive(Parser, Debug)]
#[clap(version)]
pub struct App {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Args, Debug)]
pub struct CalibrateArgs {
    /// flanking regions to omit from analysis (due to sequencing edge affects)
    #[arg(long, default_value_t = 500)]
    flank: i32,

    #[arg(short, long, default_value_t = 5678)]
    seed: u64,

    /// target fold-coverage
    #[arg(short, long, default_value_t = 40)]
    fold_coverage: i64,

    /// Size of sliding window when matching sample data coverage
    #[arg(short, long, default_value_t = 100)]
    window_size: u64,

    /// Only consider reads in the sample regions with a mapQ greater than this
    #[arg(short = 'q', long = "min-MQ", default_value_t = 10)]
    min_mapq: u8,

    /// Automatically index the output file
    #[arg(long = "write-index", default_value_t = false)]
    write_index: bool,

    /// Exclude uncalibrated (i.e., sample) reads from the output
    #[arg(short = 'x', long, default_value_t = false)]
    exclude_uncalibrated_reads: bool,

    /// Regions in the reference genome corresponding to the sequins, the name
    /// of each region must match those in the sequin BED file.
    #[arg(short = 'S', long = "sample-bed")]
    sample_bed: Option<String>,

    /// BED file specifying regions in which alignment coverage is calibrated.
    #[arg(short, long)]
    bed: String,

    /// File to write summary report (CSV). Will contain uncalibrated, target
    /// and calibrated mean coverage of each sequin region. Only created when
    /// `--sample-bed` is provided.
    #[arg(long)]
    summary_report: Option<String>,

    /// Change to experimental sample profile matching - unsuitable for
    /// production workflows.
    #[arg(long, default_value_t = false)]
    experimental: bool,

    /// Write output to file (default standard output)
    #[arg(short, long)]
    output: Option<String>,

    /// Reference sequence FASTA file. Used when input is CRAM format.
    #[arg(short = 'T', long = "reference")]
    reference: Option<String>,

    /// Write output as CRAM (requires --reference)
    #[arg(short = 'C', long = "cram", default_value_t = false)]
    cram: bool,

    path: String,
}

#[derive(Args, Debug)]
pub struct BedcovArgs {
    /// mapping quality threshold
    #[arg(short = 'Q', long = "min-MQ", default_value_t = 0)]
    min_mapq: u8,

    /// Number of bases to omit from the start and end of each region.
    #[arg(short = 'f', long = "flank", default_value_t = 0)]
    flank: i32,

    /// max per-file depth; avoids excessive memory usage. Passing zero sets it
    /// to the highest possible value, effectively removing the depth limit.
    #[arg(short = 'd', long = "max-depth", default_value_t = 8_000)]
    max_depth: u32,

    /// Reference sequence FASTA file. Used when input is CRAM format.
    #[arg(short = 'T', long = "reference")]
    reference: Option<String>,

    /// List of coverage thresholds to include in the report. The report will include the perectage
    /// of bases in each region with coverage greater than or equal to each threshold.
    #[arg(short, long, value_delimiter = ',')]
    thresholds: Option<Vec<u32>>,

    bed_path: String,
    bam_path: String,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Basic calibration of sequins
    #[command(arg_required_else_help = true)]
    Calibrate(CalibrateArgs),
    /// read depth per BED region
    Bedcov(BedcovArgs),
}

fn main() -> Result<()> {
    let args = App::parse();
    match args.command {
        Commands::Calibrate(args) => calibrate::calibrate(args)?,
        Commands::Bedcov(args) => bedcov::bedcov(args)?,
    };
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calibrate_command_parsing() {
        let args = App::parse_from([
            "app",
            "calibrate",
            "--flank",
            "1000",
            "--seed",
            "12345",
            "--fold-coverage",
            "60",
            "--window-size",
            "150",
            "--min-MQ",
            "20",
            "--write-index",
            "--exclude-uncalibrated-reads",
            "--sample-bed",
            "sample.bed",
            "--bed",
            "regions.bed",
            "--output",
            "output.txt",
            "path/to/data",
        ]);

        match args.command {
            Commands::Calibrate(calibrate_args) => {
                assert_eq!(calibrate_args.flank, 1000);
                assert_eq!(calibrate_args.seed, 12345);
                assert_eq!(calibrate_args.fold_coverage, 60);
                assert_eq!(calibrate_args.window_size, 150);
                assert_eq!(calibrate_args.min_mapq, 20);
                assert!(calibrate_args.write_index);
                assert!(calibrate_args.exclude_uncalibrated_reads);
                assert_eq!(calibrate_args.sample_bed.unwrap(), "sample.bed");
                assert_eq!(calibrate_args.bed, "regions.bed");
                assert_eq!(calibrate_args.output.unwrap(), "output.txt");
                assert_eq!(calibrate_args.path, "path/to/data");
            }
            _ => panic!("Expected Calibrate command"),
        }
    }

    #[test]
    fn test_calibrate_short_command_parsing() {
        let args = App::parse_from([
            "app",
            "calibrate",
            "--flank",
            "1000",
            "--seed",
            "12345",
            "--fold-coverage",
            "60",
            "--window-size",
            "150",
            "-q",
            "20",
            "--write-index",
            "-x",
            "-S",
            "sample.bed",
            "--bed",
            "regions.bed",
            "--output",
            "output.txt",
            "path/to/data",
        ]);

        match args.command {
            Commands::Calibrate(calibrate_args) => {
                assert_eq!(calibrate_args.flank, 1000);
                assert_eq!(calibrate_args.seed, 12345);
                assert_eq!(calibrate_args.fold_coverage, 60);
                assert_eq!(calibrate_args.window_size, 150);
                assert_eq!(calibrate_args.min_mapq, 20);
                assert!(calibrate_args.write_index);
                assert!(calibrate_args.exclude_uncalibrated_reads);
                assert_eq!(calibrate_args.sample_bed.unwrap(), "sample.bed");
                assert_eq!(calibrate_args.bed, "regions.bed");
                assert_eq!(calibrate_args.output.unwrap(), "output.txt");
                assert_eq!(calibrate_args.path, "path/to/data");
            }
            _ => panic!("Expected Calibrate command"),
        }
    }
    #[test]
    fn test_bedcov_command_parsing() {
        let args = App::parse_from([
            "app",
            "bedcov",
            "--min-MQ",
            "15",
            "--flank",
            "250",
            "regions.bed",
            "data.bam",
        ]);

        match args.command {
            Commands::Bedcov(bedcov_args) => {
                assert_eq!(bedcov_args.min_mapq, 15);
                assert_eq!(bedcov_args.flank, 250);
                assert_eq!(bedcov_args.bed_path, "regions.bed");
                assert_eq!(bedcov_args.bam_path, "data.bam");
            }
            _ => panic!("Expected Bedcov command"),
        }
    }

    #[test]
    fn test_bedcov_short_command_parsing() {
        let args = App::parse_from([
            "app",
            "bedcov",
            "-Q",
            "15",
            "-f",
            "250",
            "regions.bed",
            "data.bam",
        ]);

        match args.command {
            Commands::Bedcov(bedcov_args) => {
                assert_eq!(bedcov_args.min_mapq, 15);
                assert_eq!(bedcov_args.flank, 250);
                assert_eq!(bedcov_args.bed_path, "regions.bed");
                assert_eq!(bedcov_args.bam_path, "data.bam");
            }
            _ => panic!("Expected Bedcov command"),
        }
    }
}
