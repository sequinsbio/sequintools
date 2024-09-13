use clap::{Args, Parser, Subcommand};
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::process::exit;

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
    flank: u64,

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
        bed_path: String,
        bam_path: String,
    },
}

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

fn main() -> Result<(), Box<dyn Error>> {
    let args = App::parse();
    let _ = match args.command {
        Commands::Calibrate(cmd_args) => calibrate::calibrate(cmd_args),
        Commands::Bedcov {
            min_mapq,
            bed_path,
            bam_path,
        } => bedcov::bedcov(bam_path, bed_path, min_mapq),
    };
    Ok(())
}

fn read_bed(filename: &String) -> Vec<Region> {
    let mut vec = Vec::new();
    if let Ok(lines) = read_lines(filename) {
        for line_maybe in lines {
            if let Ok(line) = line_maybe {
                let mut contig = "";
                let mut beg = 0;
                let mut end = 0;
                let mut name = "";
                for (i, bit) in line.split_whitespace().enumerate() {
                    if i == 0 {
                        contig = bit;
                    }
                    if i == 1 {
                        beg = bit.parse::<u64>().expect(&format!(
                            "column 1 is not a integer: {}: line={}",
                            bit, line
                        ));
                    }
                    if i == 2 {
                        end = bit.parse::<u64>().expect(&format!(
                            "column 2 is not an integer: {}: line={}",
                            bit, line
                        ));
                    }
                    if i == 3 {
                        name = bit;
                    }
                }
                let r = Region {
                    contig: String::from(contig),
                    beg,
                    end,
                    name: String::from(name),
                };
                vec.push(r);
            }
        }
    }
    vec
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    match File::open(&filename) {
        Ok(file) => Ok(io::BufReader::new(file).lines()),
        Err(err) => {
            // Would be better to return an error than exiting
            println!("unable to read {}: {}", filename.as_ref().display(), err);
            exit(1);
        }
    }
}
