use anyhow::{bail, Result};
use clap::{Args, Parser, Subcommand};
use rust_htslib::bam;
use sequintools::bam::{BamReader, BamWriter, HtslibBamReader, HtslibBamWriter};
use sequintools::calibration::{self, CalibrationMode};
use sequintools::region;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(version = env!("GIT_VERSION"))]
pub struct App {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Args, Debug)]
pub struct CalibrateArgs {
    /// flanking regions to omit from analysis (due to sequencing edge affects)
    #[arg(long, default_value_t = 500)]
    flank: u64,

    #[arg(short, long, default_value_t = 5678)]
    seed: u64,

    /// target fold-coverage
    #[arg(short, long, default_value_t = 40)]
    fold_coverage: u64,

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
    sample_bed: Option<PathBuf>,

    /// BED file specifying regions in which alignment coverage is calibrated.
    #[arg(short, long)]
    bed: PathBuf,

    /// File to write summary report (CSV). Will contain uncalibrated, target
    /// and calibrated mean coverage of each sequin region. Only created when
    /// `--sample-bed` is provided.
    #[arg(long)]
    summary_report: Option<PathBuf>,

    /// Change to experimental sample profile matching - unsuitable for
    /// production workflows.
    #[arg(long, default_value_t = false)]
    experimental: bool,

    /// Write output to file (default standard output)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Reference sequence FASTA file. Used when input is CRAM format.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Write output as CRAM (requires --reference)
    #[arg(short = 'C', long = "cram", default_value_t = false)]
    cram: bool,

    path: PathBuf,
}

#[derive(Args, Debug)]
pub struct BedcovArgs {
    /// mapping quality threshold
    #[arg(short = 'Q', long = "min-MQ", default_value_t = 0)]
    min_mapq: u8,

    /// Number of bases to omit from the start and end of each region.
    #[arg(short = 'f', long = "flank", default_value_t = 0)]
    flank: u64,

    /// Reference sequence FASTA file. Used when input is CRAM format.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// List of coverage thresholds to include in the report. The report will include the
    /// percentage of bases in each region with coverage greater than or equal to each threshold.
    #[arg(short, long, value_delimiter = ',')]
    thresholds: Option<Vec<u32>>,

    bed_path: PathBuf,
    bam_path: PathBuf,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Basic calibration of sequins
    #[command(arg_required_else_help = true)]
    Calibrate(CalibrateArgs),
    /// read depth per BED region
    Bedcov(BedcovArgs),
}

impl From<BedcovArgs> for sequintools::coverage::BedcovArgs {
    fn from(args: BedcovArgs) -> Self {
        sequintools::coverage::BedcovArgs {
            min_mapq: args.min_mapq,
            flank: args.flank,
            reference: args.reference,
            thresholds: args.thresholds,
            bed_path: args.bed_path,
            bam_path: args.bam_path,
        }
    }
}

fn main() -> Result<()> {
    let args = App::parse();
    match args.command {
        Commands::Calibrate(args) => run_calibrate(&args)?,
        Commands::Bedcov(args) => sequintools::coverage::run(&args.into())?,
    };
    Ok(())
}

fn run_calibrate(args: &CalibrateArgs) -> Result<()> {
    if args.cram && args.reference.is_none() {
        bail!("--cram output requires --reference to be supplied");
    }

    let mut reader = HtslibBamReader::from_path(&args.path)?;
    let ncpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    reader.set_threads(ncpus)?;
    if let Some(reference) = args.reference.as_ref() {
        reader.set_reference(reference)?;
    }
    let hdr = bam::Header::from_template(reader.header());
    let format = if args.cram {
        bam::Format::Cram
    } else {
        bam::Format::Bam
    };

    let mut writer = if let Some(output) = &args.output {
        HtslibBamWriter::from_path(output, &hdr, format)?
    } else {
        HtslibBamWriter::from_stdout(&hdr, format)?
    };
    writer.set_threads(ncpus)?;
    if let Some(reference) = args.reference.as_ref() {
        writer.set_reference(reference)?;
    }

    let target_regions = region::load_from_bed(&mut BufReader::new(File::open(&args.bed)?))?;

    let sample_regions = if let Some(sample_bed) = &args.sample_bed {
        let regions = region::load_from_bed(&mut BufReader::new(File::open(sample_bed)?))?;
        Some(regions)
    } else {
        None
    };

    // Determine the calibration mode based on the provided arguments
    let mode = if args.experimental {
        if let Some(sample_regions) = &sample_regions {
            CalibrationMode::SampleProfile {
                sample_regions,
                flank: args.flank,
                window_size: args.window_size,
                min_mapq: args.min_mapq,
                seed: args.seed,
            }
        } else {
            return Err(anyhow::anyhow!(
                "SampleProfile mode requires sample regions. Please provide a sample BED file."
            ));
        }
    } else if let Some(sample_regions) = &sample_regions {
        CalibrationMode::SampleMeanCoverage {
            sample_regions,
            seed: args.seed,
        }
    } else {
        CalibrationMode::FixedCoverage {
            fold_coverage: args.fold_coverage,
            seed: args.seed,
        }
    };

    calibration::calibrate(&mut reader, &mut writer, &target_regions, mode)?;
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
                assert_eq!(
                    calibrate_args.sample_bed.unwrap(),
                    PathBuf::from("sample.bed")
                );
                assert_eq!(calibrate_args.bed, PathBuf::from("regions.bed"));
                assert_eq!(calibrate_args.output.unwrap(), PathBuf::from("output.txt"));
                assert_eq!(calibrate_args.path, PathBuf::from("path/to/data"));
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
                assert_eq!(
                    calibrate_args.sample_bed.unwrap(),
                    PathBuf::from("sample.bed")
                );
                assert_eq!(calibrate_args.bed, PathBuf::from("regions.bed"));
                assert_eq!(calibrate_args.output.unwrap(), PathBuf::from("output.txt"));
                assert_eq!(calibrate_args.path, PathBuf::from("path/to/data"));
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
                assert_eq!(bedcov_args.bed_path, PathBuf::from("regions.bed"));
                assert_eq!(bedcov_args.bam_path, PathBuf::from("data.bam"));
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
                assert_eq!(bedcov_args.bed_path, PathBuf::from("regions.bed"));
                assert_eq!(bedcov_args.bam_path, PathBuf::from("data.bam"));
            }
            _ => panic!("Expected Bedcov command"),
        }
    }

    #[test]
    fn test_bedcovarg_from() {
        let input = BedcovArgs {
            min_mapq: 0,
            flank: 500,
            reference: None,
            thresholds: None,
            bed_path: PathBuf::from("my.bed"),
            bam_path: PathBuf::from("my.bam"),
        };
        let expected = sequintools::coverage::BedcovArgs {
            min_mapq: 0,
            flank: 500,
            reference: None,
            thresholds: None,
            bed_path: PathBuf::from("my.bed"),
            bam_path: PathBuf::from("my.bam"),
        };
        assert_eq!(sequintools::coverage::BedcovArgs::from(input), expected);
    }
}
