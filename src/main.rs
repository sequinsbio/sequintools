use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use std::fmt;
use std::io::Read;

mod bedcov;
mod calibrate;

#[derive(Parser, Debug)]
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

    /// Write output to file (default standard output)
    #[arg(short, long)]
    output: Option<String>,

    path: String,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Basic calibration of sequins
    #[command(arg_required_else_help = true)]
    Calibrate(CalibrateArgs),
    /// read depth per BED region
    Bedcov {
        /// mapping quality threshold
        #[arg(short = 'Q', long = "min-MQ", default_value_t = 0)]
        min_mapq: u8,

        #[arg(short = 'f', long = "flank", default_value_t = 0)]
        flank: i32,

        bed_path: String,
        bam_path: String,
    },
}

#[derive(Debug, PartialEq, Eq)]
pub struct Region {
    contig: String,
    beg: u64,
    end: u64,
    name: String,
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.beg, self.end)
    }
}

fn main() -> Result<()> {
    let args = App::parse();
    let _ = match args.command {
        Commands::Calibrate(cmd_args) => calibrate::calibrate(cmd_args)?,
        Commands::Bedcov {
            min_mapq,
            flank,
            bed_path,
            bam_path,
        } => bedcov::bedcov(bam_path, bed_path, min_mapq, flank)?,
    };
    Ok(())
}

fn read_bed<R: Read>(reader: &mut R) -> Result<Vec<Region>> {
    let mut result = Vec::new();
    let mut contents = String::new();
    reader.read_to_string(&mut contents)?;
    for (i, line) in contents.lines().enumerate() {
        let bits: Vec<&str> = line.split_whitespace().collect();
        if bits.len() < 4 {
            bail!(
                "Incorrect number of columns detected, expected >= 4 found {} (line = {})",
                bits.len(),
                i + 1,
            );
        }
        let contig = bits[0];
        let beg: u64 = bits[1].parse().with_context(|| {
            format!(
                "Beg column is not an integer: is {} (line = {})",
                bits[1],
                i + 1
            )
        })?;
        let end: u64 = bits[2].parse().with_context(|| {
            format!(
                "End column is not an integer: is {} (line = {})",
                bits[2],
                i + 1
            )
        })?;
        let name = bits[3];
        result.push(Region {
            contig: String::from(contig),
            beg,
            end,
            name: String::from(name),
        });
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use std::io::{self, Read};

    #[test]
    fn test_read_bed() {
        let data = b"chr1\t1\t10\treg1";
        let mut cursor = Cursor::new(data);
        let result = read_bed(&mut cursor).unwrap();
        assert_eq!(
            result,
            vec![Region {
                contig: "chr1".to_owned(),
                beg: 1,
                end: 10,
                name: "reg1".to_owned()
            }]
        );
    }

    #[test]
    fn test_read_bed_multi() {
        let mut cursor = Cursor::new(b"chr1\t1\t10\treg1\nchr2\t2\t20\treg2");
        let result = read_bed(&mut cursor).unwrap();
        assert_eq!(
            result,
            vec![
                Region {
                    contig: "chr1".to_owned(),
                    beg: 1,
                    end: 10,
                    name: "reg1".to_owned()
                },
                Region {
                    contig: "chr2".to_owned(),
                    beg: 2,
                    end: 20,
                    name: "reg2".to_owned()
                }
            ]
        );
    }

    #[test]
    #[should_panic]
    fn test_read_bed_bad_start() {
        let mut cursor = Cursor::new(b"chr1\txxx\t10\treg1");
        read_bed(&mut cursor).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_read_bed_bad_end() {
        let mut cursor = Cursor::new(b"chr1\t1\txxx\treg1");
        read_bed(&mut cursor).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_read_bam_no_name() {
        let mut cursor = Cursor::new(b"chr1\t1\t10");
        read_bed(&mut cursor).unwrap();
    }

    struct ErrorReader;
    impl Read for ErrorReader {
        fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::Other, "bad read"))
        }
    }
    #[test]
    fn test_read_bed_bad_reader() {
        let mut reader = ErrorReader;
        let result = read_bed(&mut reader);
        assert!(result.is_err());
    }
}
