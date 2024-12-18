use crate::calibrate::mean_depth;
use crate::region::{load_from_bed, Region};
use crate::BedcovArgs;
use anyhow::Result;
use rust_htslib::bam;
use std::fs::File;
use std::io::{self, Write};

/// Generates a coverage report for genomic regions specified in a BED file.
///
/// It should calculate the mean depth via `mean_depth` function in calibrate.
///
/// This function calculates and writes coverage statistics for each region in a BED file,
/// including minimum, maximum, mean depth, standard deviation, and coefficient of variation.
///
/// # Arguments
///
/// * `bam_path` - Path to the BAM file containing aligned reads
/// * `bed_path` - Path to the BED file containing regions of interest
/// * `min_mapq` - Minimum mapping quality threshold for reads to be considered
/// * `flank` - Number of base pairs to extend regions on both sides
/// * `dest` - A writer implementing the `Write` trait where the report will be output
///
/// # Returns
///
/// Returns `Result<()>` which is:
/// * `Ok(())` if the report was successfully generated and output to the writer
/// * `Err(e)` if there was an error reading the files or writing the report
///
/// # Format
///
/// The output is a CSV file with the following columns:
/// * chrom - chromosome name
/// * beg - start position
/// * end - end position
/// * name - region name
/// * len - region length
/// * min - minimum depth
/// * max - maximum depth
/// * mean - mean depth
/// * std - standard deviation
/// * cv - coefficient of variation
fn bedcov_report<W: Write>(
    bam_path: &str,
    bed_path: &str,
    min_mapq: u8,
    flank: i32,
    dest: W,
) -> Result<()> {
    let mut bam = bam::IndexedReader::from_path(bam_path)?;
    let regions: Vec<Region> = load_from_bed(&mut io::BufReader::new(File::open(bed_path)?))?;

    let mut wtr = csv::Writer::from_writer(dest);
    wtr.write_record([
        "chrom", "beg", "end", "name", "len", "min", "max", "mean", "std", "cv",
    ])?;
    for region in &regions {
        let depth_result = mean_depth(&mut bam, region, flank, min_mapq)?;
        wtr.write_record([
            &region.contig,
            &region.beg.to_string(),
            &region.end.to_string(),
            &region.name,
            &depth_result.len().to_string(),
            &depth_result.min().unwrap_or(0).to_string(),
            &depth_result.max().unwrap_or(0).to_string(),
            &depth_result.mean().unwrap_or(0.0).to_string(),
            &depth_result.std().unwrap_or(0.0).to_string(),
            &depth_result.cv().unwrap_or(0.0).to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

/// Convenience function to generate a coverage report and write it to stdout
///
/// This is a wrapper around `bedcov_report` that writes the output to standard output.
///
/// # Arguments
///
/// * `bam_path` - Path to the BAM file containing aligned reads
/// * `bed_path` - Path to the BED file containing regions of interest
/// * `min_mapq` - Minimum mapping quality threshold for reads to be considered
/// * `flank` - Number of base pairs to extend regions on both sides
///
/// # Returns
///
/// Returns `Result<()>` which is:
/// * `Ok(())` if the report was successfully generated and written to stdout
/// * `Err(e)` if there was an error reading the files or writing the report
pub fn bedcov(args: BedcovArgs) -> Result<()> {
    bedcov_report(
        &args.bam_path,
        &args.bed_path,
        args.min_mapq,
        args.flank,
        io::stdout(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_path(key: &str) -> String {
        let file = match key {
            "bam" => "sim_R.bam",
            "bed" => "test.bed",
            _ => panic!("Unknown test file key: {}", key),
        };
        format!("{}/testdata/{}", env!("CARGO_MANIFEST_DIR"), file)
    }

    const EXPECTED_REPORT: &str = "\
    chrom,beg,end,name,len,min,max,mean,std,cv\n\
    chr1,99,199,region1,100,31,39,33.95,2.1650634,0.06377211\n\
    chr2,99,199,reg2,100,29,41,35.25,3.3087008,0.09386385";

    #[test]
    fn test_bedcov_report() {
        let mut buffer = Vec::new();
        let bam_path = &test_path("bam");
        let bed_path = &test_path("bed");
        let result = bedcov_report(bam_path, bed_path, 0, 0, &mut buffer);
        assert!(result.is_ok());

        // Compare the report with the expected
        assert_eq!(
            String::from_utf8(buffer).unwrap().trim(),
            EXPECTED_REPORT.trim()
        );
    }

    #[test]
    fn test_bedcov() {
        // BedcovArgs for testing
        let args = BedcovArgs {
            bam_path: test_path("bam"),
            bed_path: test_path("bed"),
            min_mapq: 0,
            flank: 0,
        };
        let result = bedcov(args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_bedcov_with_wrong_bam_path() {
        let args = BedcovArgs {
            bam_path: "wrong_path.bam".to_owned(),
            bed_path: test_path("bed"),
            min_mapq: 0,
            flank: 0,
        };
        let result = bedcov(args);
        assert!(result.is_err());
    }

    #[test]
    fn test_bedcov_with_invalid_bed_path() {
        let args = BedcovArgs {
            bam_path: test_path("bam"),
            bed_path: "wrong_path.bam".to_owned(),
            min_mapq: 0,
            flank: 0,
        };
        let result = bedcov(args);
        assert!(result.is_err());
    }
}
