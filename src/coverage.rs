use crate::bam::{BamReader, HtslibBamReader};
use crate::errors::{Error, Result};
use crate::region;
use crate::region::Region;
use rayon::prelude::*;
use rust_htslib::bam::record::Cigar;
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug)]
struct RegionCoverage {
    region: Region,
    coverage: Vec<u32>,
}

impl RegionCoverage {
    /// Create a new RegionCoverage instance.
    fn new(contig: &str, start: u64, end: u64, name: &str, coverage: Vec<u32>) -> Self {
        Self {
            region: Region::new(contig, start, end, name),
            coverage,
        }
    }

    /// Get the minimum coverage value.
    fn min(&self) -> Option<&u32> {
        self.coverage.iter().min()
    }

    /// Get the maximum coverage value.
    fn max(&self) -> Option<&u32> {
        self.coverage.iter().max()
    }

    /// Calculate the mean of the coverage.
    fn mean(&self) -> Option<f32> {
        let total = self.coverage.iter().sum::<u32>() as f32;
        let n = self.coverage.len() as f32;
        if n == 0.0 {
            return None;
        }
        Some(total / n)
    }

    /// Calculate the standard deviation of the coverage.
    fn std(&self) -> Option<f32> {
        match (self.mean(), self.coverage.len()) {
            (Some(mu), n) if n > 0 => {
                let variance = self
                    .coverage
                    .iter()
                    .map(|v| {
                        let diff = mu - (*v as f32);
                        diff * diff
                    })
                    .sum::<f32>()
                    / n as f32;
                Some(variance.sqrt())
            }
            _ => None,
        }
    }

    /// Calculate the coefficient of variation (CV) of the coverage.
    fn cv(&self) -> Option<f32> {
        match (self.std(), self.mean()) {
            (Some(sd), Some(mu)) if mu > 0.0 => Some(sd / mu),
            _ => None,
        }
    }

    /// Calculate the percentage of bases above a certain coverage threshold.
    fn percent_above_threshold(&self, threshold: u32) -> Option<f64> {
        let n = self.coverage.len();
        if n == 0 {
            return None;
        }
        let x = self.coverage.iter().filter(|&&x| x >= threshold).count();
        Some(x as f64 / n as f64)
    }
}

fn coverage_for_region<T: BamReader>(
    bam_reader: &mut T,
    region: &Region,
    min_mapq: u8,
    flank: u64,
    reference: Option<&PathBuf>,
) -> Result<RegionCoverage> {
    let ncpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    bam_reader.set_threads(ncpus)?;
    if let Some(reference) = reference {
        bam_reader.set_reference(reference)?;
    }
    let beg = region.beg.checked_add(flank).ok_or_else(|| Error::Bedcov {
        msg: format!(
            "Region start + flank overflowed for region: {}:{}-{}",
            region.contig, region.beg, region.end
        ),
    })?;
    let end = region.end.checked_sub(flank).ok_or_else(|| Error::Bedcov {
        msg: format!(
            "Region end - flank underflowed for region: {}:{}-{}",
            region.contig, region.beg, region.end
        ),
    })?;
    if beg >= end {
        return Err(Error::Bedcov {
            msg: format!(
                "Region start >= end after applying flank for region: {}:{}-{}",
                region.contig, region.beg, region.end
            ),
        });
    }

    let mut coverage = vec![0u32; (end - beg) as usize];

    let tid = bam_reader
        .header()
        .tid(region.contig.as_bytes())
        .ok_or_else(|| Error::Bedcov {
            msg: format!("Chromosome {} not found in BAM header", region.contig),
        })?;
    bam_reader.fetch((tid, beg, end))?;

    for result in bam_reader.records() {
        let record = result?;
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        if record.mapq() < min_mapq {
            continue;
        }
        let read_start = record.pos();
        let mut ref_pos = read_start;

        for &cigar_op in record.cigar().iter() {
            match cigar_op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    for i in 0..len {
                        let pos = ref_pos + i as i64;
                        if pos >= beg as i64 && pos < end as i64 {
                            let idx = (pos - beg as i64) as usize;
                            coverage[idx] += 1;
                        }
                    }
                    ref_pos += len as i64
                }
                Cigar::Del(len) | Cigar::RefSkip(len) => ref_pos += len as i64,
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                    // These do not consume reference positions
                }
            }
        }
    }
    Ok(RegionCoverage::new(
        region.contig.as_str(),
        region.beg,
        region.end,
        region.name.as_str(),
        coverage,
    ))
}

fn calculate_coverage(
    bam_path: &PathBuf,
    regions: &[Region],
    reference: Option<&PathBuf>,
    min_mapq: u8,
    flank: u64,
) -> Result<Vec<RegionCoverage>> {
    let result = regions
        .par_iter()
        .map(|region| {
            // rust_htslib::bam::IndexedReader is not Send + Sync (thread
            // safe). Each thread needs its own copy (I think).
            let mut bam_reader = HtslibBamReader::from_path(bam_path)?;
            let coverage =
                coverage_for_region(&mut bam_reader, region, min_mapq, flank, reference)?;
            Ok(coverage)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(result)
}

/// Write coverage results to CSV format.
///
/// This function writes the coverage results to a CSV file, including
/// optional thresholds.
///
/// # Example
/// ```ignore
/// let coverages = vec![
///     RegionCoverage::new("chr1", 100, 200, "region1", vec![1, 2, 3]),
///     RegionCoverage::new("chr1", 200, 300, "region2", vec![4, 5, 6]),
/// ];
/// let thresholds = Some(vec![2, 4]);
/// write_csv(&coverages, thresholds, std::io::stdout())?;
/// ```
fn write_csv<W: Write>(
    coverages: &[RegionCoverage],
    thresholds: Option<Vec<u32>>,
    mut dest: W,
) -> Result<()> {
    let mut columns: Vec<String> = [
        "name", "chrom", "beg", "end", "min", "max", "mean", "std", "cv",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    if let Some(thresholds) = &thresholds {
        for threshold in thresholds {
            columns.push(format!("pct_gt_{threshold}"));
        }
    }
    writeln!(dest, "{}", columns.join(","))?;
    for coverage in coverages.iter() {
        let min = coverage.min().ok_or_else(|| Error::Bedcov {
            msg: "unable to calculate minimum".to_string(),
        })?;
        let max = coverage.max().ok_or_else(|| Error::Bedcov {
            msg: "unable to calculate maximum".to_string(),
        })?;
        let mean = coverage.mean().ok_or_else(|| Error::Bedcov {
            msg: "unable to calculate mean".to_string(),
        })?;
        let std = coverage.std().ok_or_else(|| Error::Bedcov {
            msg: "unable to calculate std".to_string(),
        })?;
        let cv = coverage.cv().ok_or_else(|| Error::Bedcov {
            msg: "unable to calculate cv".to_string(),
        })?;
        let mut row = format!(
            "{},{},{},{},{min},{max},{mean:.2},{std:.2},{cv:.4}",
            coverage.region.name, coverage.region.contig, coverage.region.beg, coverage.region.end,
        );
        if let Some(thresholds) = &thresholds {
            for thresh in thresholds {
                let pct = coverage.percent_above_threshold(*thresh);
                match pct {
                    Some(value) => row.push_str(&format!(",{:.2}", value)),
                    None => {
                        return Err(Error::Bedcov {
                            msg: format!("unable to calculate pct_gt_{}", thresh),
                        })
                    }
                }
            }
        }
        writeln!(dest, "{row}")?;
    }
    Ok(())
}

#[derive(Debug, PartialEq, Eq)]
pub struct BedcovArgs {
    pub bam_path: PathBuf,
    pub bed_path: PathBuf,
    pub reference: Option<PathBuf>,
    pub min_mapq: u8,
    pub flank: u64,
    pub thresholds: Option<Vec<u32>>,
}

pub fn run(args: &BedcovArgs) -> Result<()> {
    let file = std::fs::File::open(&args.bed_path)?;
    let mut reader = std::io::BufReader::new(file);
    let regions = region::load_from_bed(&mut reader)?;
    let coverages = calculate_coverage(
        &args.bam_path,
        &regions,
        args.reference.as_ref(),
        args.min_mapq,
        args.flank,
    )?;
    write_csv(&coverages, args.thresholds.clone(), std::io::stdout())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::MockBamReader;
    use rust_htslib::bam::record::CigarString;
    use rust_htslib::bam::Record;

    #[test]
    fn test_write_csv() {
        let coverages = vec![
            RegionCoverage::new("chr1", 100, 200, "region1", vec![1, 2, 3]),
            RegionCoverage::new("chr1", 200, 300, "region2", vec![4, 5, 6]),
        ];
        let mut output = Vec::new();
        let result = write_csv(&coverages, None, &mut output);
        assert!(result.is_ok());

        let expected = "\
name,chrom,beg,end,min,max,mean,std,cv
region1,chr1,100,200,1,3,2.00,0.82,0.4082
region2,chr1,200,300,4,6,5.00,0.82,0.1633";
        assert_eq!(String::from_utf8(output).unwrap().trim(), expected.trim());
    }

    #[test]
    fn test_write_csv_with_thresholds() {
        let coverages = vec![
            RegionCoverage::new("chr1", 100, 200, "region1", vec![1, 2, 3]),
            RegionCoverage::new("chr1", 200, 300, "region2", vec![4, 5, 6]),
        ];
        let thresholds = Some(vec![2, 4]);
        let mut output = Vec::new();
        let result = write_csv(&coverages, thresholds, &mut output);
        assert!(result.is_ok());

        let expected = "\
name,chrom,beg,end,min,max,mean,std,cv,pct_gt_2,pct_gt_4
region1,chr1,100,200,1,3,2.00,0.82,0.4082,0.67,0.00
region2,chr1,200,300,4,6,5.00,0.82,0.1633,1.00,1.00";
        assert_eq!(String::from_utf8(output).unwrap().trim(), expected.trim());
    }

    #[test]
    fn test_write_csv_no_coverage() {
        let coverages = vec![RegionCoverage::new("chr1", 100, 200, "region1", vec![])];
        let mut output = Vec::new();
        let result = write_csv(&coverages, None, &mut output);
        assert!(result.is_err());
        let err = result.unwrap_err();
        match err {
            Error::Bedcov { msg } => {
                assert_eq!(msg, "unable to calculate minimum");
            }
            _ => panic!("Expected Bedcov error"),
        }
    }

    #[test]
    fn test_coverage_for_region_empty() {
        let records = vec![];
        let mut mock = MockBamReader::new(records, None);
        let region = Region::new("chrQ_mirror", 100, 200, "test_region");
        let result = coverage_for_region(&mut mock, &region, 0, 0, None);
        assert!(result.is_ok());
        let coverage = result.unwrap();
        let max = coverage.max().expect("should have max");
        assert_eq!(*max, 0);
    }

    #[test]
    fn test_coverage_for_region() {
        let mut record = Record::new();
        record.set_pos(100);
        record.set_cigar(Some(&CigarString(vec![Cigar::Match(100)])));
        record.set_mapq(60);
        record.set_flags(99); // PAIRED,PROPER_PAIR,MREVERSE,READ1
        let mut mock = MockBamReader::new(vec![record], None);
        let region = Region::new("chrQ_mirror", 100, 200, "test_region");
        let result = coverage_for_region(&mut mock, &region, 0, 0, None);
        assert!(result.is_ok());
        let coverage = result.unwrap();
        let max = coverage.max().expect("should have max");
        assert_eq!(*max, 1);
        assert_eq!(coverage.mean().unwrap(), 1.0);
    }

    #[test]
    fn test_coverage_for_region_missing_chrom() {
        let mut record = Record::new();
        record.set_pos(100);
        record.set_cigar(Some(&CigarString(vec![Cigar::Match(100)])));
        record.set_mapq(60);
        record.set_flags(99); // PAIRED,PROPER_PAIR,MREVERSE,READ1
        let mut mock = MockBamReader::new(vec![record], None);
        let region = Region::new("chr2", 100, 200, "test_region");
        let result = coverage_for_region(&mut mock, &region, 0, 0, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_coverage_for_region_skip_reads() {
        let mut record1 = Record::new();
        record1.set_pos(100);
        record1.set_cigar(Some(&CigarString(vec![Cigar::Match(100)])));
        record1.set_mapq(60);
        record1.set_flags(99); // PAIRED,PROPER_PAIR,MREVERSE,READ1

        let mut record2 = Record::new();
        record2.set_pos(100);
        record2.set_flags(68); // READ1,UNMAP

        let mut record3 = Record::new();
        record3.set_pos(100);
        record3.set_flags(99); // PAIRED,PROPER_PAIR,MREVERSE,READ1
        record3.set_mapq(0);

        let mut mock = MockBamReader::new(vec![record1, record2, record3], None);
        let region = Region::new("chrQ_mirror", 100, 200, "test_region");
        let min_mapq = 20;
        let result = coverage_for_region(&mut mock, &region, min_mapq, 0, None);
        assert!(result.is_ok());
        let coverage = result.unwrap();
        let max = coverage.max().expect("should have max");
        assert_eq!(*max, 1);
        assert_eq!(coverage.mean().unwrap(), 1.0);
    }

    #[test]
    fn test_calculate_coverage() {
        let bam_path = PathBuf::from("testdata/calibrated.bam");
        let region = Region::new("chrQ_mirror", 100, 200, "test_region");
        let result = calculate_coverage(&bam_path, &[region], None, 0, 0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_run() {
        let args = BedcovArgs {
            bam_path: PathBuf::from("testdata/calibrated.bam"),
            bed_path: PathBuf::from("testdata/calibrated.bed"),
            reference: None,
            min_mapq: 0,
            flank: 0,
            thresholds: None,
        };
        let result = run(&args);
        assert!(result.is_ok(), "Failed to run bedcov: {result:?}");
    }

    #[test]
    fn test_coverage_for_region_flank_overflow() {
        let records = vec![];
        let mut mock = MockBamReader::new(records, None);
        let region = Region::new("chrQ_mirror", 100, 200, "test_region");
        let flank = u64::MAX;
        let result = coverage_for_region(&mut mock, &region, 0, flank, None);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err
            .to_string()
            .contains("Region start + flank overflowed for region"));
    }

    #[test]
    fn test_coverage_for_region_flank_underflow() {
        let records = vec![];
        let mut mock = MockBamReader::new(records, None);
        let region = Region::new("chrQ_mirror", u64::MIN, 200, "test_region");
        let flank = u64::MAX - 10;
        let result = coverage_for_region(&mut mock, &region, 0, flank, None);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err
            .to_string()
            .contains("Region end - flank underflowed for region"));
    }

    #[test]
    fn test_coverage_for_region_flank_bigger_than_region() {
        let records = vec![];
        let mut mock = MockBamReader::new(records, None);
        let region = Region::new("chrQ_mirror", 100, 200, "test_region");
        let flank = 150;
        let result = coverage_for_region(&mut mock, &region, 0, flank, None);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err
            .to_string()
            .contains("Region start >= end after applying flank for region"));
    }
}
