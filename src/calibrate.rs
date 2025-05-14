//! Calibration module for sequencing data analysis.
//!
//! This module provides functionality for calibrating sequencing data based on coverage analysis.
//! It supports two main calibration methods:
//! - Standard coverage calibration: Applies a mean target coverage to all sequin regions
//! - Sample coverage calibration: Calibrates based on sample-specific coverage patterns
//!
//! # Main Functions
//!
//! - [`calibrate`]: Main entry point for calibration operations
//! - [`calibrate_by_standard_coverage`]: Performs standard coverage-based calibration
//! - [`calibrate_by_sample_coverage`]: Performs sample-specific coverage calibration
//! - [`mean_depth`]: Calculates mean depth for a given region
//!
//! # Types
//!
//! - [`DepthResult`]: Represents depth calculation results including histogram and statistics
//! - [`CalibrateError`]: Custom error type for calibration operations
//!
//! # Examples
//!
//! ```no_run
//! use sequintools::CalibrateArgs;
//! use sequintools::calibrate;
//!
//! let args = CalibrateArgs {
//!     // ... configure arguments
//! };
//!
//! match calibrate::calibrate(args) {
//!     Ok(_) => println!("Calibration successful"),
//!     Err(e) => eprintln!("Calibration failed: {}", e),
//! }
//! ```
//!
//! # Features
//!
//! - BAM file processing and indexing
//! - Region-specific depth calculation
//! - Coverage histogram generation
//! - Statistical analysis of coverage data
//! - Support for paired-end reads
//! - Configurable mapping quality filters
//!
//! # Technical Details
//!
//! The module handles both standard and sample-based calibration methods:
//! - Standard calibration applies a uniform coverage target
//! - Sample-based calibration uses sliding windows to match coverage patterns
//! - Both methods preserve read pairs and handle various read flags
//!
//! # Note
//!
//! Ensure input BAM files are properly indexed before calibration.
//! The module assumes 0-based coordinates in BED files.
use crate::region::{self, Region};
use crate::CalibrateArgs;
use anyhow::{bail, Context, Result};
use rand::seq::IteratorRandom;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;
use rust_htslib::bam::{IndexedReader, Writer};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io;
use std::path::Path;
use std::process::exit;
use std::vec::Vec;

// Defines the maximum depth for pileup operations. htslib defines this as
// 2_147_483_647 which is the maximum value of i32. The rust_htslib function
// this is passed to requires an u32 hence the cast.
const PILEUP_MAX_DEPTH: u32 = i32::MAX as u32;

/// Calibrates sequin coverage based on the provided arguments.
///
/// This function serves as the main entry point for calibration operations. It determines the
/// appropriate calibration method (standard or sample-specific) based on the presence of a sample
/// BED file in the arguments and delegates the calibration process to the corresponding function.
///
/// # Arguments
///
/// * `args` - [`CalibrateArgs`] containing:
///   - `bed`: Path to BED file with regions to calibrate
///   - `path`: Path to input BAM file
///   - `output`: Optional output BAM path (uses stdout if None)
///   - `sample_bed`: Optional path to sample BED file with regions to use for calibration
///   - `fold_coverage`: Target coverage fold
///   - `flank`: Number of bases to exclude from region edges
///   - `seed`: Random seed for reproducible sampling
///   - `write_index`: Whether to write BAM index
///
/// # Returns
///
/// Returns `Result<()>` which is:
/// - `Ok(())` if calibration completed successfully
/// - `Err(CalibrateError)` if an error occurs during calibration
///
/// # Errors
///
/// Will return error if:
/// - Cannot open input BAM file
/// - Cannot create output BAM file
/// - Region processing or depth calculation fails
///
/// # Example
///
/// ```no_run
/// use sequintools::CalibrateArgs;
/// use sequintools::calibrate;
///
/// let args = CalibrateArgs {
///     bed: "regions.bed".to_string(),
///     path: "input.bam".to_string(),
///     output: Some("output.bam".to_string()),
///     sample_bed: None,
///     fold_coverage: 30,
///     // ... other fields
/// };
///
/// match calibrate(args) {
///     Ok(_) => println!("Calibration successful"),
///     Err(e) => eprintln!("Calibration failed: {}", e),
/// }
/// ```
pub fn calibrate(args: CalibrateArgs) -> Result<()> {
    if !args.experimental {
        calibrate_by_standard_coverage(args)?;
    } else {
        calibrate_by_sample_coverage(args)?;
    }
    Ok(())
}

/// Calibrates sequin coverage by applying a mean target coverage to all sequin
/// regions.
///
/// The target coverage is either obtained from the mean coverage of the
/// corresponding region in the sample data (if `args.sample_bed` is not None),
/// or using a fixed coverage (taken from `args.fold_coverage`).
///
/// This function performs standard coverage calibration by:
/// 1. Processing each region in the input BED file
/// 2. Calculating mean depth for each region
/// 3. Randomly sampling reads to achieve the target coverage
/// 4. Preserving read pairs in the output
///
/// # Arguments
///
/// * `args` - [`CalibrateArgs`] containing:
///   - `bed`: Path to BED file with regions to calibrate
///   - `path`: Path to input BAM file
///   - `output`: Optional output BAM path (uses stdout if None)
///   - `fold_coverage`: Target coverage fold
///   - `flank`: Number of bases to exclude from region edges
///   - `seed`: Random seed for reproducible sampling
///   - `write_index`: Whether to write BAM index
///
/// # Returns
///
/// Returns `Result<()>` which is:
/// - `Ok(())` if calibration completed successfully
/// - `Err(CalibrateError)` if sample bed file is provided or other errors occur
///
/// # Errors
///
/// Will return error if:
/// - Sample BED file is provided (incompatible with standard calibration)
/// - Cannot open input BAM file
/// - Cannot create output BAM file
/// - Region processing or depth calculation fails
///
/// # Example
///
/// ```no_run
/// use sequintools::CalibrateArgs;
///
/// let args = CalibrateArgs {
///     bed: "regions.bed".to_string(),
///     path: "input.bam".to_string(),
///     output: Some("output.bam".to_string()),
///     sample_bed: None,
///     fold_coverage: 30,
///     // ... other fields
/// };
///
/// calibrate_by_standard_coverage(args)?;
/// ```
pub fn calibrate_by_standard_coverage(args: CalibrateArgs) -> Result<()> {
    // Check summary report file doesn't already exist, so that we don't spend time calibrating
    // only to fail later.
    if let Some(path) = args.summary_report.as_ref() {
        if Path::new(&path).exists() {
            bail!("The summary report file '{}' already exists", path);
        }
    };

    let regions = region::load_from_bed(&mut io::BufReader::new(File::open(&args.bed)?))?;

    let sample_regions = args
        .sample_bed
        .as_ref()
        .map(File::open)
        .transpose()?
        .map(|mut file| region::load_from_bed(&mut file))
        .transpose()?
        .map(|regions| {
            regions
                .into_iter()
                .map(|region| (region.name.clone(), region))
                .collect::<HashMap<_, _>>()
        });

    let mut bam = match bam::IndexedReader::from_path(&args.path) {
        Ok(r) => r,
        Err(err) => {
            eprintln!("unable to open input BAM: {}", err);
            exit(1);
        }
    };
    let header = bam::Header::from_template(bam.header());
    let mut out = match &args.output {
        Some(path) => bam::Writer::from_path(path, &header, bam::Format::Bam)?,
        None => bam::Writer::from_stdout(&header, bam::Format::Bam)?,
    };

    // This assumes that the decoy chromosome is the last one in the BAM file.
    // If there are multiple decoy chromosomes, the behavior is not well
    // defined.

    if !args.exclude_uncalibrated_reads {
        copy_uncalibrated_contigs(&mut bam, &mut out, &regions)?;
    }

    let mut calibration_results =
        calibrate_regions_by_fixed_coverage(&mut bam, &mut out, &regions, &sample_regions, &args)?;

    if !args.exclude_uncalibrated_reads {
        copy_unmapped_reads(&mut bam, &mut out)?;
    }

    drop(out);

    if args.write_index {
        if let Some(path) = &args.output {
            bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
        }
    }

    if args.summary_report.is_some() {
        write_summary_report(&mut calibration_results, &args)?;
    }

    Ok(())
}

/// Create a CSV summary of the calibration results.
fn write_summary_report(
    calibration_results: &mut [CalibrationResult],
    args: &CalibrateArgs,
) -> Result<()> {
    let path = args.summary_report.as_ref().unwrap();
    let file = File::create(path)?;
    let mut writer = csv::Writer::from_writer(file);
    writer.write_record([
        "name",
        "chrom",
        "start",
        "end",
        "uncalibrated_coverage",
        "target_coverage",
        "calibrated_coverage",
    ])?;

    for result in calibration_results {
        let region = &result.region;
        // How can we compute the calibrated coverage if no file is written (i.e., it's being
        // piped)
        if let Some(path) = &args.output {
            let mut bam = IndexedReader::from_path(path)?;
            let flank = if args.sample_bed.is_some() {
                0
            } else {
                args.flank
            };
            result.calibrated_coverage =
                mean_depth(&mut bam, region, flank, args.min_mapq, PILEUP_MAX_DEPTH)?.mean;
        }

        writer.write_record(&[
            region.name.clone(),
            region.contig.clone(),
            region.beg.to_string(),
            region.end.to_string(),
            result.uncalibrated_coverage.to_string(),
            result.target_coverage.to_string(),
            result.calibrated_coverage.to_string(),
        ])?;
    }
    Ok(())
}

fn copy_uncalibrated_contigs(
    bam: &mut bam::IndexedReader,
    out: &mut bam::Writer,
    regions: &[Region],
) -> Result<()> {
    let hdr = bam.header().clone();
    let calibrated_contigs = regions
        .iter()
        .map(|r| r.contig.clone())
        .collect::<HashSet<String>>();
    let calibrated_tids = (0..hdr.target_count())
        .filter(|tid| {
            calibrated_contigs.contains(&String::from_utf8_lossy(hdr.tid2name(*tid)).into_owned())
        })
        .map(|tid| tid as i32)
        .collect::<Vec<i32>>();

    eprintln!("Copying uncalibrated reads");
    for tid in 0..hdr.target_count() {
        if calibrated_tids.contains(&(tid as i32)) {
            continue;
        }
        bam.fetch(tid)?;
        for result in bam.records() {
            let record = result?;
            if !calibrated_tids.contains(&record.tid()) {
                out.write(&record)?;
            }
        }
    }
    Ok(())
}

struct CalibrationResult {
    region: Region,
    uncalibrated_coverage: f64,
    target_coverage: f64,
    calibrated_coverage: f64,
}

fn calibrate_regions_by_fixed_coverage(
    bam: &mut IndexedReader,
    out: &mut Writer,
    regions: &[Region],
    sample_regions: &Option<HashMap<String, Region>>,
    args: &CalibrateArgs,
) -> Result<Vec<CalibrationResult>> {
    let mut rng = Pcg32::seed_from_u64(args.seed);
    let mut keep_names = HashMap::new();

    // If we're using the sample mean coverage, we don't want to exclude the flanks from any of the
    // calculations.
    let flank = match args.sample_bed {
        Some(_) => 0,
        None => args.flank,
    };
    let mut results = Vec::new();
    for region in regions {
        let target_coverage = sample_regions
            .as_ref()
            .and_then(|hash| hash.get(&region.name))
            .and_then(|region| mean_depth(bam, region, flank, args.min_mapq, PILEUP_MAX_DEPTH).ok())
            .map(|result| result.mean)
            .unwrap_or(args.fold_coverage as f64);

        let depth = mean_depth(bam, region, flank, args.min_mapq, PILEUP_MAX_DEPTH)?.mean;
        let threshold = target_coverage / depth;

        eprintln!(
            "Calibrating {} ({}:{}-{}) mean_coverage={} target_coverage={}",
            region.name, region.contig, region.beg, region.end, depth, target_coverage
        );

        bam.fetch((&region.contig, region.beg, region.end))?;
        for r in bam.records() {
            let record = r?;
            if subsample(&record, &mut keep_names, threshold, &mut rng) {
                out.write(&record)?;
            }
        }
        let result = CalibrationResult {
            region: region.clone(),
            uncalibrated_coverage: depth,
            target_coverage,
            calibrated_coverage: 0.0,
        };
        results.push(result);
    }
    Ok(results)
}

fn subsample(
    record: &bam::Record,
    hash: &mut HashMap<String, bool>,
    threshold: f64,
    rng: &mut Pcg32,
) -> bool {
    if record.is_duplicate() {
        return false;
    }
    let qname = String::from_utf8(record.qname().to_vec()).unwrap();
    match hash.get(&qname) {
        Some(_) => return true,
        None => {
            if record.pos() > record.mpos() {
                // We've already considered the mate, and didn't keep it. We
                // can't keep this read as it will create a singleton.
                return false;
            }
        }
    };
    let rand = rng.random::<f64>();
    if rand <= threshold {
        hash.insert(qname, true);
        true
    } else {
        false
    }
}

fn copy_unmapped_reads(bam: &mut bam::IndexedReader, out: &mut bam::Writer) -> Result<()> {
    bam.fetch("*")?;
    for result in bam.records() {
        let record = result?;
        out.write(&record)?;
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

/// Return the mean depth of a region from a BAM file. `flank` bases are removed
/// from the start and end of the region prior to depth calculation. The region
/// must use 0-based coordinates.
///
/// # Arguments
///
/// * `bam` - A mutable reference to an indexed BAM reader.
/// * `region` - The region for which to calculate mean depth.
/// * `flank` - Number of bases to exclude from the start and end of the region.
/// * `min_mapq` - Minimum mapping quality for reads to be considered.
/// * `max_depth` - Maximum depth for pile. Setting this to 0 uses the maximum possible value, effectively removing the depth limit.
///
/// # Returns
///
/// Returns `Result<DepthResult>` which is:
/// - `Ok(DepthResult)` if depth calculation completed successfully
/// - `Err(anyhow::Error)` if an error occurs during depth calculation
///
/// # Errors
///
/// Will return error if:
/// - Cannot fetch region from BAM file
/// - Pileup operation fails
///
/// # Example
///
/// ```no_run
/// use sequintools::Region;
/// use sequintools::calibrate::mean_depth;
/// use rust_htslib::bam;
///
/// let mut bam = bam::IndexedReader::from_path("input.bam").unwrap();
/// let region = Region {
///     contig: "chr1".to_string(),
///     beg: 100,
///     end: 200,
///     name: "region1".to_string(),
/// };
///
/// let depth_result = mean_depth(&mut bam, &region, 10, 20, 8_000).unwrap();
/// println!("Mean depth: {}", depth_result.mean);
/// ```
pub fn mean_depth(
    bam: &mut bam::IndexedReader,
    region: &Region,
    flank: i32,
    min_mapq: u8,
    max_depth: u32,
) -> Result<DepthResult> {
    // If max_depth is 0, set it to the maximum possible value, effectively
    // removing the depth limit. Note, internally htslib uses i32 for depth so
    // we cannot exceed i32::MAX despite it being passed as u32.
    let max_depth = if max_depth == 0 {
        i32::MAX as u32
    } else {
        max_depth
    };

    let beg: u32 = region.beg as u32 + flank as u32;
    let end: u32 = region.end as u32 - flank as u32;
    bam.fetch((region.contig.as_str(), beg, end))?;
    let mut total = 0;
    let len = end - beg;

    let include_flags: u16 = 0;
    // 0xf04 3844 UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
    let exclude_flags: u16 = 3844;

    let mut xs: HashMap<u32, u32> = std::collections::HashMap::new();

    let mut pileups = bam.pileup();
    pileups.set_max_depth(max_depth);

    for p in pileups {
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

/// Calibrates sequin coverage based on sample-specific coverage patterns.
///
/// This function performs sample-specific coverage calibration by:
/// 1. Loading sample regions from the provided BED file
/// 2. Iterating over each contig in the BAM file
/// 3. For each contig, if it is a calibrated contig, it calibrates the regions using sample data
/// 4. If the contig is not calibrated and `exclude_uncalibrated_reads` is false, it processes all reads in the contig
/// 5. Writes the calibrated reads to the output BAM file
///
/// # Arguments
///
/// * `args` - [`CalibrateArgs`] containing:
///   - `bed`: Path to BED file with regions to calibrate
///   - `path`: Path to input BAM file
///   - `output`: Optional output BAM path (uses stdout if None)
///   - `sample_bed`: Path to sample BED file with regions to use for calibration
///   - `flank`: Number of bases to exclude from region edges
///   - `seed`: Random seed for reproducible sampling
///   - `window_size`: Size of the sliding window for coverage calculation
///   - `min_mapq`: Minimum mapping quality for reads to be considered
///   - `write_index`: Whether to write BAM index
///   - `exclude_uncalibrated_reads`: Whether to exclude reads from uncalibrated contigs
///
/// # Returns
///
/// Returns `Result<()>` which is:
/// - `Ok(())` if calibration completed successfully
/// - `Err(CalibrateError)` if sample bed file is not provided or other errors occur
///
/// # Errors
///
/// Will return error if:
/// - Sample BED file is not provided (required for sample-specific calibration)
/// - Cannot open input BAM file
/// - Cannot create output BAM file
/// - Region processing or depth calculation fails
///
/// # Example
///
/// ```no_run
/// use sequintools::CalibrateArgs;
///
/// let args = CalibrateArgs {
///     bed: "regions.bed".to_string(),
///     path: "input.bam".to_string(),
///     output: Some("output.bam".to_string()),
///     sample_bed: Some("sample.bed".to_string()),
///     // ... other fields
/// };
///
/// calibrate_by_sample_coverage(args)?;
/// ```
fn calibrate_by_sample_coverage(args: CalibrateArgs) -> Result<()> {
    // Should NOT enter this function if sample_bed in args is not provided.
    if args.sample_bed.is_none() {
        return Err(CalibrateError::new(
            "Cannot use sample coverage calibration when sample_bed file is not provided.",
        )
        .into());
    }

    let cal_contigs = calibrated_contigs(&args.bed)?;

    let regions = region::load_from_bed(&mut io::BufReader::new(File::open(&args.bed)?))?;
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

        let mut sample_regions_map = HashMap::new();
        let sample_regions = region::load_from_bed(&mut io::BufReader::new(File::open(
            args.sample_bed.as_ref().unwrap(),
        )?))?;
        for region in &sample_regions {
            sample_regions_map.insert(&region.name, region);
        }

        for i in 0..hdr.target_count() {
            let contig = target_names[i as usize];
            if cal_contigs.contains(contig) {
                calibrate_regions(&mut bam, &mut out, &regions, &sample_regions_map, &args);
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

/// Calibrates the reads in a set of regions by randomly sampling reads from the sample data.
///
/// The sample data is assumed to be a BED file with regions that have been sequenced. The sample
/// data is used to determine the coverage in each window of the region. The coverage in the region
/// is then adjusted to match the coverage in the sample data by randomly sampling reads that start
/// in each window of the region. The number of reads to sample is determined by the coverage in the
/// sample data. Both ends of paired-end reads are kept.
///
/// # Arguments
///
/// * `bam` - A mutable reference to an indexed BAM reader.
/// * `out` - A mutable reference to a BAM writer for the output.
/// * `regions` - A vector of regions to calibrate.
/// * `sample_regions_map` - A hashmap mapping region names to sample regions.
/// * `args` - Calibration arguments containing various parameters.
///
/// # Errors
///
/// This function does not return errors directly but may panic if there are issues with BAM file
/// operations or if required regions are not found in the sample regions map.
fn calibrate_regions(
    bam: &mut bam::IndexedReader,
    out: &mut bam::Writer,
    regions: &Vec<Region>,
    sample_regions_map: &HashMap<&String, &Region>,
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
                Err(_) => (0..region_records.len())
                    .map(|x| x.try_into().unwrap())
                    .collect(),
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
    let (beg, end) = (region.beg as i64, region.end as i64);
    let mut n: i64 = 0;
    bam.fetch((&region.contig, beg, end)).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        if record.pos() >= beg && record.pos() <= end && record.mapq() >= min_mapq {
            n += 1;
        }
    }
    n
}

/// Calculates the number of read starts in sliding windows across a specified region.
///
/// This function iterates over the specified region in the BAM file, dividing it into sliding windows
/// of a given size. For each window, it counts the number of reads that start within that window
/// and returns these counts as a vector.
///
/// # Arguments
///
/// * `bam` - A mutable reference to an indexed BAM reader.
/// * `region` - The region for which to calculate read starts.
/// * `flank` - Number of bases to exclude from the start and end of the region.
/// * `window_size` - The size of each sliding window.
/// * `min_mapq` - Minimum mapping quality for reads to be considered.
///
/// # Returns
///
/// A vector of integers where each element represents the number of read starts in a corresponding window.
pub fn window_starts(
    bam: &mut bam::IndexedReader,
    region: &Region,
    flank: u64,
    window_size: u64,
    min_mapq: u8,
) -> Vec<i64> {
    let (contig, name) = (&region.contig, &region.name);
    let mut starts = Vec::new();

    let region_beg = region.beg + flank;
    let region_end = region.end - flank - window_size;

    for beg in (region_beg..region_end).step_by(window_size as usize) {
        let end = beg + window_size - 1;
        let n = starts_in(
            bam,
            &Region {
                contig: contig.to_owned(),
                beg,
                end,
                name: name.to_owned(),
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
            details: msg.to_owned(),
        }
    }
}

impl fmt::Display for CalibrateError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.details)
    }
}

// Implement the Error trait to make it a proper error type
impl Error for CalibrateError {}

fn choose_from(size: i64, n: i64, seed: u64) -> Result<Vec<i64>, CalibrateError> {
    let mut rng = Pcg32::seed_from_u64(seed);
    let xs = (0..size).choose_multiple(&mut rng, n as usize);
    if xs.len() as i64 == n {
        Ok(xs)
    } else {
        let msg = format!("sample size is too big, wanted {} from {}", n, size);
        Err(CalibrateError::new(&msg))
    }
}

// Return the unique set of contig names in a BED file
fn calibrated_contigs(path: &str) -> Result<HashSet<String>> {
    let mut contigs = HashSet::new();
    let regions = region::load_from_bed(&mut io::BufReader::new(File::open(path)?))?;
    for region in &regions {
        contigs.insert(region.contig.to_owned());
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
    use tempfile::NamedTempFile;

    const TEST_BAM_PATH: &str = "testdata/sim_R.bam";
    const TEST_BED_PATH: &str = "testdata/region_test.bed";

    fn mock_calibrate_args(with_sample: bool, need_output: bool) -> CalibrateArgs {
        CalibrateArgs {
            bed: TEST_BED_PATH.to_owned(),
            path: TEST_BAM_PATH.to_owned(),
            output: if need_output {
                // Create temporary file for output, it ensures a unique file name.
                Some(
                    NamedTempFile::new()
                        .unwrap()
                        .path()
                        .to_string_lossy()
                        .into_owned(),
                )
            } else {
                None
            },
            sample_bed: if with_sample {
                Some(TEST_BED_PATH.to_owned())
            } else {
                None
            },
            flank: 0,
            seed: 42,
            fold_coverage: 1,
            window_size: 10,
            min_mapq: 0,
            write_index: false,
            exclude_uncalibrated_reads: false,
            experimental: false,
            summary_report: None,
        }
    }

    #[test]
    fn test_write_summary_report() {
        let mut args = mock_calibrate_args(false, false);
        args.summary_report = Some(
            NamedTempFile::new()
                .unwrap()
                .path()
                .to_string_lossy()
                .into_owned(),
        );
        let mut calibration_results = vec![
            CalibrationResult {
                region: Region {
                    contig: "chrQ".to_string(),
                    beg: 1,
                    end: 2,
                    name: "region1".to_string(),
                },
                uncalibrated_coverage: 1.0,
                target_coverage: 2.0,
                calibrated_coverage: 3.0,
            },
            CalibrationResult {
                region: Region {
                    contig: "chrQ".to_string(),
                    beg: 10,
                    end: 20,
                    name: "region2".to_string(),
                },
                uncalibrated_coverage: 4.0,
                target_coverage: 5.0,
                calibrated_coverage: 6.0,
            },
        ];
        let result = write_summary_report(&mut calibration_results, &args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_write_summary_report_with_output() {
        let mut args = mock_calibrate_args(false, false);
        args.output = Some("testdata/calibrated.bam".to_string());
        args.summary_report = Some(
            NamedTempFile::new()
                .unwrap()
                .path()
                .to_string_lossy()
                .into_owned(),
        );
        let mut calibration_results = vec![
            CalibrationResult {
                region: Region {
                    contig: "chrQ_mirror".to_string(),
                    beg: 200,
                    end: 4342,
                    name: "SG_000000038".to_string(),
                },
                uncalibrated_coverage: 1.0,
                target_coverage: 2.0,
                calibrated_coverage: 3.0,
            },
            CalibrationResult {
                region: Region {
                    contig: "chrQ_mirror".to_string(),
                    beg: 4542,
                    end: 8684,
                    name: "SG_000000039".to_string(),
                },
                uncalibrated_coverage: 4.0,
                target_coverage: 5.0,
                calibrated_coverage: 6.0,
            },
        ];
        let result = write_summary_report(&mut calibration_results, &args);
        assert!(result.is_ok());

        let mut rdr = csv::Reader::from_path(args.summary_report.as_ref().unwrap()).unwrap();
        let header = rdr.headers().unwrap();
        assert_eq!(
            header,
            &csv::StringRecord::from(vec![
                "name",
                "chrom",
                "start",
                "end",
                "uncalibrated_coverage",
                "target_coverage",
                "calibrated_coverage"
            ])
        );
        let mut records = rdr.records();
        let r1 = records.next().unwrap().unwrap();
        assert_eq!(
            r1,
            csv::StringRecord::from(vec![
                "SG_000000038",
                "chrQ_mirror",
                "200",
                "4342",
                "1",
                "2",
                "0"
            ])
        );
        let r2 = records.next().unwrap().unwrap();
        assert_eq!(
            r2,
            csv::StringRecord::from(vec![
                "SG_000000039",
                "chrQ_mirror",
                "4542",
                "8684",
                "4",
                "5",
                "0"
            ])
        );
        assert!(records.next().is_none());
    }

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
            PILEUP_MAX_DEPTH,
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
        let result = mean_depth(&mut bam, &region, 0, 0, PILEUP_MAX_DEPTH)
            .unwrap()
            .mean;
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
        let result = mean_depth(&mut bam, &region, 0, 10, PILEUP_MAX_DEPTH)
            .unwrap()
            .mean;
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_mean_depth_increased_max_depth() {
        let mut bam = bam::IndexedReader::from_path("testdata/bedcov_max_depth.bam").unwrap();
        let region = Region {
            contig: "chrQ_mirror".to_owned(),
            beg: 500,
            end: 800,
            name: "test_region".to_owned(),
        };
        let result = mean_depth(&mut bam, &region, 0, 0, 100_000_000).unwrap();
        assert_eq!(result.max().unwrap(), 11444);
    }

    #[test]
    fn test_mean_depth_decreased_max_depth() {
        let mut bam = bam::IndexedReader::from_path("testdata/bedcov_max_depth.bam").unwrap();
        let region = Region {
            contig: "chrQ_mirror".to_owned(),
            beg: 500,
            end: 800,
            name: "test_region".to_owned(),
        };
        let result = mean_depth(&mut bam, &region, 0, 0, 1_000).unwrap();
        assert_eq!(result.max().unwrap(), 1121);
    }

    #[test]
    fn test_mean_depth_max_depth_0() {
        let mut bam = bam::IndexedReader::from_path("testdata/bedcov_max_depth.bam").unwrap();
        let region = Region {
            contig: "chrQ_mirror".to_owned(),
            beg: 500,
            end: 800,
            name: "test_region".to_owned(),
        };
        let result = mean_depth(&mut bam, &region, 0, 0, 0).unwrap();
        assert_eq!(result.max().unwrap(), 11444);
    }

    #[test]
    fn test_calibrated_error() {
        let err = CalibrateError::new("test error");
        assert_eq!(err.to_string(), "test error");
    }

    #[test]
    fn test_choose_from() {
        // Test successful case
        let result = choose_from(10, 5, 42);
        assert!(result.is_ok());
        let values = result.unwrap();
        assert_eq!(values.len(), 5);
        assert!(values.iter().all(|&x| (0..10).contains(&x)));

        // Test requesting more items than available
        let result = choose_from(5, 10, 42);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "sample size is too big, wanted 10 from 5"
        );

        // Test empty source range
        let result = choose_from(0, 1, 42);
        assert!(result.is_err());

        // Test requesting zero items
        let result = choose_from(5, 0, 42);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);

        // Test deterministic results with same seed
        let result1 = choose_from(10, 3, 42).unwrap();
        let result2 = choose_from(10, 3, 42).unwrap();
        assert_eq!(result1, result2);
    }

    #[test]
    fn test_calibrate_regions() {
        // Setup test data
        let mut bam = bam::IndexedReader::from_path(TEST_BAM_PATH).unwrap();
        let header = bam::Header::from_template(bam.header());

        let regions =
            region::load_from_bed(&mut io::BufReader::new(File::open(TEST_BED_PATH).unwrap()))
                .unwrap();
        let mut sample_regions_map = HashMap::new();
        for region in &regions {
            sample_regions_map.insert(&region.name, region);
        }
        let args = mock_calibrate_args(true, true);
        let out_path = args.output.clone().unwrap();
        let mut temp_output = bam::Writer::from_path(&out_path, &header, bam::Format::Bam).unwrap();

        // Run calibration
        calibrate_regions(
            &mut bam,
            &mut temp_output,
            &regions,
            &sample_regions_map,
            &args,
        );

        // Flush and close output
        drop(temp_output);

        // Verify output file exists and has content
        let metadata = std::fs::metadata(&out_path).unwrap();
        assert!(metadata.len() > 0, "Output BAM file should not be empty");

        // TODO: check the output BAM
    }

    #[test]
    fn test_calibrate_by_sample_coverage() {
        // Throw error if args don't provide sample
        let args_without_sample = mock_calibrate_args(false, false);
        let result = calibrate_by_sample_coverage(args_without_sample);
        assert!(result.is_err());

        // Working correctly with expected args
        let args_expected = mock_calibrate_args(true, true);
        let out_path = args_expected.output.clone().unwrap();

        let result = calibrate_by_sample_coverage(args_expected);
        assert!(result.is_ok());

        // Verify output file exists and has content
        let metadata = std::fs::metadata(&out_path).unwrap();
        assert!(metadata.len() > 0, "Output BAM file should not be empty");
    }

    #[test]
    fn test_calibrate_by_standard_coverage() {
        // Working correctly with expected args
        let args_expected = mock_calibrate_args(false, true);
        let out_path = args_expected.output.clone().unwrap();

        let result = calibrate_by_standard_coverage(args_expected);
        assert!(result.is_ok());

        // Verify output file exists and has content
        let metadata = std::fs::metadata(&out_path).unwrap();
        assert!(metadata.len() > 0, "Output BAM file should not be empty");

        // TODO: check the output BAM
    }

    #[test]
    fn test_calibrate_by_standard_coverage_with_sample() {
        let args_expected = mock_calibrate_args(true, true);
        let result = calibrate_by_standard_coverage(args_expected);
        assert!(result.is_ok());
    }

    #[test]
    fn test_subsample() {
        let mut rng = Pcg32::seed_from_u64(42);
        let mut hash = HashMap::new();
        let mut record = bam::Record::new();
        record.set_qname(b"read1");
        record.set_pos(100);
        record.set_mpos(200);

        let result = subsample(&record, &mut hash, 1.0, &mut rng);
        assert!(result);

        // Test with key added by first subsample
        let result = subsample(&record, &mut hash, 0.0, &mut rng);
        assert!(result);
    }

    #[test]
    fn test_subsample_reject_read() {
        let mut rng = Pcg32::seed_from_u64(42);
        let mut hash = HashMap::new();
        let mut record = bam::Record::new();
        record.set_qname(b"read1");
        record.set_pos(100);
        record.set_mpos(200);

        // Test with a different threshold
        let result = subsample(&record, &mut hash, 0.0, &mut rng);
        assert!(!result);
    }

    #[test]
    fn test_subsample_mate_before_read() {
        let mut rng = Pcg32::seed_from_u64(42);
        let mut hash = HashMap::new();
        let mut record = bam::Record::new();
        record.set_qname(b"read1");
        record.set_pos(200);
        record.set_mpos(100);
        let result = subsample(&record, &mut hash, 1.0, &mut rng);
        assert!(!result);
    }

    #[test]
    fn test_subsample_read_in_hash() {
        let mut rng = Pcg32::seed_from_u64(42);

        let mut hash = HashMap::new();
        hash.insert("read1".to_string(), true);

        let mut record = bam::Record::new();
        record.set_qname(b"read1");
        record.set_pos(100);
        record.set_mpos(200);

        let result = subsample(&record, &mut hash, 1.0, &mut rng);
        assert!(result);
    }

    #[test]
    fn test_calibrate() {
        let args = mock_calibrate_args(false, true);
        let result = calibrate(args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_calibrate_experimental() {
        let mut args = mock_calibrate_args(true, true);
        args.experimental = true;
        let result = calibrate(args);
        assert!(result.is_ok(), "{:?}", result);
    }

    #[test]
    fn test_calibrate_experimental_missing_sample_bed() {
        let mut args = mock_calibrate_args(false, true);
        args.experimental = true;
        let result = calibrate(args);
        assert!(result.is_err(), "{:?}", result);
    }
}
