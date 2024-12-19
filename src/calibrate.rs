use crate::region::{load_from_bed, Region};
use crate::CalibrateArgs;
use anyhow::{Context, Result};
use rand::seq::IteratorRandom;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io;
use std::process::exit;
use std::vec::Vec;

pub fn calibrate(args: CalibrateArgs) -> Result<()> {
    let _ = match args.sample_bed {
        Some(_) => calibrate_by_sample_coverage(args),
        None => calibrate_by_standard_coverage(args),
    };
    Ok(())
}

// calibrate sequin coverage by applying a mean target coverage to all sequin
// regions.
pub fn calibrate_by_standard_coverage(args: CalibrateArgs) -> Result<()> {
    let f = File::open(&args.bed)?;
    let mut reader = io::BufReader::new(f);
    let regions = load_from_bed(&mut reader)?;

    let mut bam = match bam::IndexedReader::from_path(&args.path) {
        Ok(r) => r,
        Err(err) => {
            eprintln!("unable to open input BAM: {}", err);
            exit(1);
        }
    };
    {
        let header = bam::Header::from_template(bam.header());
        let mut out = match &args.output {
            Some(path) => bam::Writer::from_path(path, &header, bam::Format::Bam)?,
            None => bam::Writer::from_stdout(&header, bam::Format::Bam)?,
        };

        let mut rng = Pcg32::seed_from_u64(args.seed);

        let mut pairs = HashMap::new();
        let coverage = args.fold_coverage as f64;

        for region in &regions {
            let depth = mean_depth(&mut bam, region, args.flank, 10)?.mean;
            println!(
                "{}:{}-{} {} {}",
                region.contig,
                region.beg,
                region.end,
                depth,
                depth * (coverage / depth)
            );

            let chrom = region.contig.as_str();
            let beg = region.beg + args.flank as u64;
            let end = region.end - args.flank as u64;

            bam.fetch((chrom, beg, end)).unwrap();

            for r in bam.records() {
                let record = r.unwrap();
                let qname = String::from_utf8(record.qname().to_vec()).unwrap();

                match pairs.get(&qname) {
                    Some(_) => out.write(&record).unwrap(),
                    _ => {
                        let n = rng.gen::<f64>();
                        let thres = coverage / depth;
                        if n <= thres {
                            // sequintools: thres = coverage / mean; if rand <= keep record
                            // need to ensure we keep both reads of a pair
                            out.write(&record).unwrap();
                            pairs.insert(qname, true);
                        }
                    }
                }
            }
        }
    }
    if args.write_index {
        if let Some(path) = &args.output {
            bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
        }
    }
    Ok(())
}

pub struct DepthResult {
    pub histogram: Vec<(u32, u32)>,
    pub mean: f64,
    len: u32,
}

impl DepthResult {
    pub fn cv(&self) -> Option<f32> {
        match (self.std(), self.mean()) {
            (Some(sig), Some(mu)) => Some(sig / mu),
            _ => None,
        }
    }

    pub fn len(&self) -> i32 {
        self.len as i32
    }

    pub fn min(&self) -> Option<i32> {
        self.histogram.iter().map(|x| x.1 as i32).min()
    }

    pub fn max(&self) -> Option<i32> {
        self.histogram.iter().map(|x| x.1 as i32).max()
    }

    pub fn mean(&self) -> Option<f32> {
        let data: Vec<i32> = self.histogram.iter().map(|x| x.1 as i32).collect();
        let sum = data.iter().sum::<i32>() as f32;
        let count = data.len();
        match count {
            positive if positive > 0 => Some(sum / count as f32),
            _ => None,
        }
    }

    pub fn std(&self) -> Option<f32> {
        let data: Vec<i32> = self.histogram.iter().map(|x| x.1 as i32).collect();
        match (self.mean(), data.len()) {
            (Some(data_mean), count) if count > 0 => {
                let variance = data
                    .iter()
                    .map(|value| {
                        let diff = data_mean - (*value as f32);
                        diff * diff
                    })
                    .sum::<f32>()
                    / count as f32;
                Some(variance.sqrt())
            }
            _ => None,
        }
    }
}

// Return the mean depth of a region from a BAM file. `flank` bases are removed
// from the start and end of the region prior to depth calculation. The region
// must use 0-based coordinates.
pub fn mean_depth(
    bam: &mut bam::IndexedReader,
    region: &Region,
    flank: i32,
    min_mapq: u8,
) -> Result<DepthResult> {
    let beg: u32 = region.beg as u32 + flank as u32;
    let end: u32 = region.end as u32 - flank as u32;
    bam.fetch((region.contig.as_str(), beg, end))?;
    let mut total = 0;
    let len = end - beg;

    let include_flags: u16 = 0;
    // 0xf04 3844 UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
    let exclude_flags: u16 = 3844;

    let mut xs: HashMap<u32, u32> = std::collections::HashMap::new();

    for p in bam.pileup() {
        let pileup = p.with_context(|| "pileup failed".to_string())?;
        if pileup.pos() < beg || pileup.pos() > end - 1 {
            continue;
        }

        let depth = pileup
            .alignments()
            .filter(|aln| {
                let record = aln.record();
                let flags = record.flags();
                (!flags) & include_flags == 0
                    && flags & exclude_flags == 0
                    && record.mapq() >= min_mapq
            })
            .count() as u32;

        xs.insert(pileup.pos(), depth);
        total += depth;
    }
    let mut hist = vec![(0, 0); len as usize];
    for (i, pos) in (beg..end).enumerate() {
        hist[i] = (pos, xs.get(&pos).copied().unwrap_or(0));
    }
    Ok(DepthResult {
        histogram: hist,
        mean: total as f64 / len as f64,
        len,
    })
}

//  What if we use a sliding window over the region and determine the depth of
//  that window in the sample data. For the sequin data we randomly sample but
//  only if the read starts in that window. We keep both ends of a pair. How do
//  we handle secondary/supplementry/duplicate reads? We can experiment with
//  different window sizes to get the best replication of the sample coverage.
fn calibrate_by_sample_coverage(args: CalibrateArgs) -> Result<()> {
    // shouldn't be calling this function if sample_bed is not provided
    let sample_bed = args.sample_bed.as_ref().unwrap();
    let cal_contigs = calibrated_contigs(&args.bed)?;

    let f = File::open(&args.bed)?;
    let mut reader = io::BufReader::new(f);
    let regions = load_from_bed(&mut reader)?;
    let mut bam = match bam::IndexedReader::from_path(&args.path) {
        Ok(r) => r,
        Err(err) => {
            eprintln!("unable to open input BAM: {}", err);
            exit(1);
        }
    };
    {
        let header = bam::Header::from_template(bam.header());
        let mut out = match &args.output {
            Some(path) => bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap(),
            None => bam::Writer::from_stdout(&header, bam::Format::Bam).unwrap(),
        };

        let hdr = bam.header().clone();

        let target_lengths: Vec<u64> = (0..hdr.target_count())
            .map(|i| hdr.target_len(i).unwrap())
            .collect();

        let target_names: Vec<&str> = hdr
            .target_names()
            .into_iter()
            .map(|name| std::str::from_utf8(name).unwrap())
            .collect();

        for i in 0..hdr.target_count() {
            let contig = target_names[i as usize];
            if cal_contigs.contains(contig) {
                let mut sample_regions_map = HashMap::new();
                let f = File::open(sample_bed)?;
                let mut reader = io::BufReader::new(f);
                let sample_regions = load_from_bed(&mut reader)?;
                for region in &sample_regions {
                    sample_regions_map.insert(&region.name, region);
                }
                calibrate_regions(&mut bam, &mut out, &regions, sample_regions_map, &args);
            } else if !args.exclude_uncalibrated_reads {
                eprintln!("Processing reads in {}", contig);
                let beg = 1;
                let end = target_lengths[i as usize];
                bam.fetch((contig, beg, end)).unwrap();
                for r in bam.records() {
                    let record = r.unwrap();
                    out.write(&record).unwrap();
                }
            }
        }
    }
    if args.write_index {
        if let Some(path) = &args.output {
            eprintln!("Writing index");
            bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
        }
    }
    Ok(())
}

fn calibrate_regions(
    bam: &mut bam::IndexedReader,
    out: &mut bam::Writer,
    regions: &Vec<Region>,
    sample_regions_map: HashMap<&String, &Region>,
    args: &CalibrateArgs,
) {
    // iterate over sequin regions
    for region in regions {
        let contig = &region.contig;
        let name = &region.name;

        eprintln!(
            "Calibrating region {} {}:{}-{}",
            name, contig, region.beg, region.end
        );
        let sample_region = sample_regions_map.get(name).unwrap();

        let sample_starts = window_starts(
            bam,
            sample_region,
            args.flank as u64,
            args.window_size,
            args.min_mapq,
        );
        let rev_sample_starts: Vec<i64> = sample_starts.into_iter().rev().collect();

        let region_beg = region.beg + args.flank as u64;
        let region_end = region.end - args.flank as u64 - args.window_size;

        // This should be all reads mapped to a region
        let records = records_that_start_in_region(bam, contig, region.beg, region.end);
        let mut keep_names = HashSet::new();

        for (i, window_beg) in (region_beg..region_end)
            .step_by(args.window_size as usize)
            .enumerate()
        {
            let window_end = window_beg + args.window_size - 1;

            // We are always keeping both reads from paired-end data, so
            // ever  ytime we select a read to keep we are actually keeping two
            // reads. This results in higher coverage in the calibrated sequin
            // regions compared to the sample region. I'm dividing by 2 as a
            // quick and dirty adjustment; it gives a sequin coverage that is
            // much closer to the sample region (it is actually slightly lower
            // than the sample coverage).
            let n_starts = rev_sample_starts[i] / 2;

            // These are the regions that *start* in the current window
            let region_records: Vec<&bam::Record> = records
                .iter()
                .filter(|record| {
                    let pos: u64 = record.pos().try_into().unwrap();
                    pos >= window_beg && pos <= window_end
                })
                .collect();

            let indexes = choose_from(region_records.len() as i64, n_starts, args.seed);
            let numbers = match indexes {
                Ok(numbers) => numbers,
                Err(_) => {
                    // This only seems to happen in the last window, why?
                    // eprintln!(
                    //     "warning: not enough records while processing {} {}:{}-{} window {}-{}: {}, using all records",
                    //     name, contig, region_beg, region_end, window_beg, window_end, err
                    // );
                    (0..region_records.len())
                        .map(|x| x.try_into().unwrap())
                        .collect()
                }
            };
            for idx in &numbers {
                let record = region_records[*idx as usize];
                let qname = String::from_utf8(record.qname().to_vec()).unwrap();
                keep_names.insert(qname);
            }
        }

        for record in &records {
            let qname = String::from_utf8(record.qname().to_vec()).unwrap(); // really? there must be a better way
            if keep_names.contains(&qname) {
                out.write(record).unwrap();
            }
        }
    }
}

fn starts_in(bam: &mut bam::IndexedReader, region: &Region, min_mapq: u8) -> i64 {
    let chrom = region.contig.as_str();
    let beg = region.beg as i64;
    let end = region.end as i64;
    let mut n: i64 = 0;
    bam.fetch((chrom, beg, end)).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        if record.pos() < beg || record.pos() > end {
            continue;
        }
        if record.mapq() < min_mapq {
            continue;
        }
        n += 1;
    }
    n
}

pub fn window_starts(
    bam: &mut bam::IndexedReader,
    region: &Region,
    flank: u64,
    window_size: u64,
    min_mapq: u8,
) -> Vec<i64> {
    let contig = &region.contig;
    let name = &region.name;
    let region_beg = region.beg + flank;
    let region_end = region.end - flank - window_size;
    let mut starts = Vec::new();

    for beg in (region_beg..region_end).step_by(window_size as usize) {
        let end = beg + window_size - 1;
        let n = starts_in(
            bam,
            &Region {
                contig: contig.to_string(),
                beg,
                end,
                name: name.to_string(),
            },
            min_mapq,
        );
        starts.push(n);
    }
    starts
}

fn records_that_start_in_region(
    bam: &mut bam::IndexedReader,
    contig: &str,
    beg: u64,
    end: u64,
) -> Vec<bam::Record> {
    let mut records = Vec::new();
    bam.fetch((contig, beg, end)).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        let pos = record.pos() as u64;
        if pos < beg || pos > end {
            continue;
        }
        records.push(record);
    }
    records
}

#[derive(Debug)]
struct CalibrateError {
    details: String,
}

impl CalibrateError {
    fn new(msg: &str) -> CalibrateError {
        CalibrateError {
            details: msg.to_string(),
        }
    }
}

impl fmt::Display for CalibrateError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.details)
    }
}

impl Error for CalibrateError {
    fn description(&self) -> &str {
        &self.details
    }
}

fn choose_from(size: i64, n: i64, seed: u64) -> Result<Vec<i64>, CalibrateError> {
    let mut rng = Pcg32::seed_from_u64(seed);
    let xs = (0..size).choose_multiple(&mut rng, n as usize);
    if xs.len() as i64 == n {
        Ok(xs)
    } else {
        let msg = format!("sample size is too big, wanted {} from {}", n, size);
        Err(CalibrateError::new(msg.as_str()))
    }
}

// Return the unique set of contig names in a BED file
fn calibrated_contigs(path: &str) -> Result<HashSet<String>> {
    let mut contigs = HashSet::new();
    let f = File::open(path)?;
    let mut reader = io::BufReader::new(f);
    let regions = load_from_bed(&mut reader)?;
    for region in &regions {
        let contig = String::from(&region.contig);
        contigs.insert(contig);
    }
    Ok(contigs)
}

#[cfg(test)]
impl DepthResult {
    pub fn create_for_test(histogram: Vec<(u32, u32)>, mean: f64, len: u32) -> Self {
        // Provide a default instance for testing
        Self {
            histogram,
            mean,
            len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_BAM_PATH: &str = "testdata/sim_R.bam";
    const TEST_BED_PATH: &str = "testdata/region_test.bed";

    #[test]
    fn test_window_starts() {
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let region = Region {
            contig: "chr1".to_owned(),
            beg: 99,
            end: 199,
            name: "region".to_owned(),
        };
        let starts = window_starts(&mut bam, &region, 0, 10, 0);
        assert_eq!(starts.len(), 9);
    }

    #[test]
    fn test_window_starts_with_mapq() {
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let region = Region {
            contig: "chr3".to_owned(),
            beg: 99,
            end: 878,
            name: "reg2".to_owned(),
        };
        let starts = window_starts(&mut bam, &region, 0, 100, 98);

        let expected = [36, 37, 35, 36, 20, 29, 32];
        for (i, cnt) in starts.iter().enumerate() {
            assert_eq!(
                cnt, &expected[i],
                "Window {}: Expected {} but Got {} ",
                i, expected[i], cnt
            );
        }
    }

    #[test]
    fn test_records_that_start_in_region() {
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let records = records_that_start_in_region(&mut bam, "chr3", 99, 199);
        assert_eq!(records.len(), 36);
    }

    #[test]
    fn test_calibrated_contigs() {
        let contigs = calibrated_contigs(TEST_BED_PATH).unwrap();
        assert!(contigs.contains("chr1"));
        assert!(contigs.contains("chr2"));
        assert_eq!(contigs.len(), 2);
    }

    #[test]
    fn test_mean_depth() {
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let result = mean_depth(
            &mut bam,
            &Region {
                contig: "chr1".to_owned(),
                beg: 99,
                end: 199,
                name: "region".to_owned(),
            },
            0,
            0,
        );
        assert_eq!(result.unwrap().mean, 33.95);
    }

    #[test]
    fn test_mean_depth_chr2() {
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let region = Region {
            contig: "chr2".to_owned(),
            beg: 99,
            end: 199,
            name: "reg2".to_owned(),
        };
        let result = mean_depth(&mut bam, &region, 0, 0).unwrap().mean;
        assert_eq!(result, 35.25);
    }

    #[test]
    fn test_mean_depth_chr2_mapq10() {
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let region = Region {
            contig: "chr2".to_owned(),
            beg: 99,
            end: 199,
            name: "reg2".to_owned(),
        };
        let result = mean_depth(&mut bam, &region, 0, 10).unwrap().mean;
        assert_eq!(result, 0.0);
    }
}
