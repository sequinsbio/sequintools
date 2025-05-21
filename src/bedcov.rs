use crate::calibrate::{mean_depth, DepthResult};
use crate::region::{load_from_bed, Region};
use crate::BedcovArgs;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io::{self, Write};

const REPORT_HEADER: [&str; 10] = [
    "chrom", "beg", "end", "name", "len", "min", "max", "mean", "std", "cv",
];

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
/// * `chrom` - Chromosome name.
/// * `beg` - Start position of the region.
/// * `end` - End position of the region.
/// * `name` - Name of the region.
/// * `len` - Length of the region.
/// * `min` - Minimum depth within the region.
/// * `max` - Maximum depth within the region.
/// * `mean` - Mean depth across the region.
/// * `std` - Standard deviation of the depth.
/// * `cv` - Coefficient of variation of the depth.
fn bedcov_report<W: Write>(
    bam_path: &str,
    bed_path: &str,
    min_mapq: u8,
    flank: i32,
    max_depth: u32,
    reference_path: Option<String>,
    dest: W,
) -> Result<()> {
    let mut bam = bam::IndexedReader::from_path(bam_path)?;
    let ncpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    bam.set_threads(ncpus)?;
    if let Some(reference) = reference_path {
        bam.set_reference(reference)?;
    }
    let regions: Vec<Region> = load_from_bed(&mut io::BufReader::new(File::open(bed_path)?))?;

    let mut wtr = csv::Writer::from_writer(dest);
    wtr.write_record(REPORT_HEADER)?;
    for region in &regions {
        let depth_result = mean_depth(&mut bam, region, flank, min_mapq, max_depth)?;
        let record = build_line_for_header(region, &depth_result)?;
        wtr.write_record(record)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Builds a line of CSV data corresponding to the given region and depth result.
///
/// This function constructs a vector of strings representing a single line of CSV output,
/// with each element in the vector corresponding to a column in the CSV file. The columns
/// are defined by the `REPORT_HEADER` constant, and the values are extracted from the
/// provided `region` and `depth_result`.
///
/// # Arguments
///
/// * `region` - A reference to a `Region` struct representing a genomic region.
/// * `depth_result` - A reference to a `DepthResult` struct containing depth statistics for the region.
///
/// # Returns
///
/// Returns a `Result` containing:
/// * `Ok(Vec<String>)` - A vector of strings representing the CSV line if successful.
/// * `Err(Error)` - An error if an unexpected header is encountered.
///
/// # Errors
///
/// This function will return an error if an unexpected header is encountered in the `REPORT_HEADER` array.
fn build_line_for_header(region: &Region, depth_result: &DepthResult) -> Result<Vec<String>> {
    let mut line = vec!["".to_string(); REPORT_HEADER.len()];
    for (i, header) in REPORT_HEADER.iter().enumerate() {
        line[i] = match *header {
            "chrom" => region.contig.clone(),
            "beg" => region.beg.to_string(),
            "end" => region.end.to_string(),
            "name" => region.name.clone(),
            "len" => depth_result.len().to_string(),
            "min" => depth_result.min().unwrap_or(0).to_string(),
            "max" => depth_result.max().unwrap_or(0).to_string(),
            "mean" => depth_result.mean().unwrap_or(0.0).to_string(),
            "std" => depth_result.std().unwrap_or(0.0).to_string(),
            "cv" => depth_result.cv().unwrap_or(0.0).to_string(),
            _ => panic!("Unexpected header: {}", header),
        };
    }
    Ok(line)
}

/// Generates a coverage report for genomic regions and outputs it to the standard output.
///
/// This is a convenient wrapper around the `bedcov_report` function. It takes the necessary
/// parameters to generate the coverage report and writes the resulting CSV directly to
/// the standard output (stdout).
///
/// # Arguments
///
/// * `args` - A `BedcovArgs` struct wrapping args from commandline, which containing:
///     * `bam_path` - Path to the BAM file with aligned reads.
///     * `bed_path` - Path to the BED file with regions of interest.
///     * `min_mapq` - Minimum mapping quality score for reads to be included.
///     * `flank` - Number of base pairs to extend each region on both sides.
///
/// # Returns
///
/// A `Result<()>` indicating the success or failure of the report generation:
///
/// * `Ok(())` - The report was successfully generated and written to stdout.
/// * `Err(e)` - An error occurred during file processing or output operations.
pub fn bedcov(args: BedcovArgs) -> Result<()> {
    bedcov_report(
        &args.bam_path,
        &args.bed_path,
        args.min_mapq,
        args.flank,
        args.max_depth,
        args.reference,
        io::stdout(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_line() {
        let region = Region {
            contig: "chr1".to_owned(),
            beg: 100,
            end: 200,
            name: "region1".to_owned(),
        };

        // Create a mock DepthResult for testing using direct struct initialization
        let depth_result = DepthResult::create_for_test(
            vec![(1, 33), (1, 31), (1, 36), (2, 35), (2, 37)],
            30.0,
            10,
        );
        let line = build_line_for_header(&region, &depth_result).unwrap();

        assert_eq!(
            line,
            vec![
                "chr1",
                "100",
                "200",
                "region1",
                "10",
                "31",
                "37",
                "34.4",
                "2.154066",
                "0.062618196",
            ]
        );
    }

    fn get_test_path(key: &str) -> String {
        let file = match key {
            "bam" => "sim_R.bam",
            "bed" => "region_test.bed",
            _ => panic!("Unknown test file key: {}", key),
        };
        format!("{}/testdata/{}", env!("CARGO_MANIFEST_DIR"), file)
    }

    const EXPECTED_REPORT: &str = "\
    chrom,beg,end,name,len,min,max,mean,std,cv\n\
    chr1,99,199,region1,100,31,39,33.95,2.1650634,0.06377211\n\
    chr2,99,199,reg2,100,29,41,35.25,3.3087008,0.09386385";

    #[test]
    fn bedcov_report_with_wtr() {
        let mut buffer = Vec::new();
        let bam_path = &get_test_path("bam");
        let bed_path = &get_test_path("bed");
        let reference_path = None;
        let result = bedcov_report(bam_path, bed_path, 0, 0, 8_000, reference_path, &mut buffer);
        assert!(result.is_ok());

        // Compare the report with the expected
        assert_eq!(
            String::from_utf8(buffer).unwrap().trim(),
            EXPECTED_REPORT.trim()
        );
    }

    #[test]
    fn bedcov_stdout() {
        // BedcovArgs for testing
        let args = BedcovArgs {
            bam_path: get_test_path("bam"),
            bed_path: get_test_path("bed"),
            min_mapq: 0,
            flank: 0,
            max_depth: 8_000,
            reference: None,
        };
        let result = bedcov(args);
        assert!(result.is_ok());
    }

    #[test]
    fn bedcov_with_wrong_bam_path() {
        let args = BedcovArgs {
            bam_path: "wrong_path.bam".to_owned(),
            bed_path: get_test_path("bed"),
            min_mapq: 0,
            flank: 0,
            max_depth: 8_000,
            reference: None,
        };
        let result = bedcov(args);
        assert!(result.is_err());
    }

    #[test]
    fn bedcov_with_invalid_bed_path() {
        let args = BedcovArgs {
            bam_path: get_test_path("bam"),
            bed_path: "wrong_path.bam".to_owned(),
            min_mapq: 0,
            flank: 0,
            max_depth: 8_000,
            reference: None,
        };
        let result = bedcov(args);
        assert!(result.is_err());
    }

    const TEST_CRAM_PATH: &str = "testdata/calibrated.cram";
    const TEST_CRAM_REF_PATH: &str = "testdata/reference.fasta";
    const TEST_CRAM_BED_PATH: &str = "testdata/cram_regions.bed";

    #[test]
    fn bedcov_stdout_cram() {
        // BedcovArgs for testing
        let args = BedcovArgs {
            bam_path: TEST_CRAM_PATH.to_string(),
            bed_path: TEST_CRAM_BED_PATH.to_string(),
            min_mapq: 0,
            flank: 0,
            max_depth: 8_000,
            reference: Some(TEST_CRAM_REF_PATH.to_string()),
        };
        let result = bedcov(args);
        assert!(result.is_ok());
    }
}
