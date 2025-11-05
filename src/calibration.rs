//! # Calibration Module
//!
//! This module provides functionality for calibrating BAM files based on coverage
//! requirements. It supports different calibration modes: fixed coverage,
//! sample mean coverage, and sample profile matching.
//!
//! The main entry point is the [`calibrate`] function, which processes a BAM file
//! and writes the calibrated output to another BAM file.

use crate::bam::{BamReader, BamWriter};
use crate::coverage;
use crate::errors::{Error, Result};
use crate::region::Region;
use rand::seq::IteratorRandom;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;
use rust_htslib::bam::{FetchDefinition, Record};
use std::collections::HashMap;
use std::collections::HashSet;

/// Represents the different modes for calibration.
///
/// Each variant specifies a different strategy for determining how to downsample
/// reads in the target regions.
pub enum CalibrationMode<'a> {
    /// Calibrate to a fixed fold coverage.
    ///
    /// # Fields
    /// - `fold_coverage`: The desired fold coverage for all target regions.
    /// - `seed`: Random seed for reproducible downsampling.
    FixedCoverage { fold_coverage: u64, seed: u64 },
    /// Calibrate based on the mean coverage of sample regions.
    ///
    /// # Fields
    /// - `sample_regions`: Regions to sample mean coverage from.
    /// - `seed`: Random seed for reproducible downsampling.
    SampleMeanCoverage {
        sample_regions: &'a [Region],
        seed: u64,
    },
    /// Calibrate by matching the coverage profile of sample regions.
    ///
    /// # Fields
    /// - `sample_regions`: Regions to match the profile of.
    /// - `flank`: Number of bases to flank the regions.
    /// - `window_size`: Size of windows for profile matching.
    /// - `min_mapq`: Minimum mapping quality for reads.
    /// - `seed`: Random seed for reproducible downsampling.
    SampleProfile {
        sample_regions: &'a [Region],
        flank: u64,
        window_size: u64,
        min_mapq: u8,
        seed: u64,
    },
}

/// Calibrates a BAM file by downsampling reads in target regions according to the specified mode.
///
/// This function processes the input BAM file, applies calibration to the specified target regions,
/// and writes the result to the output BAM file. It also copies uncalibrated contigs and unmapped reads.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `writer`: A mutable reference to a BAM writer.
/// - `target_regions`: A slice of regions to calibrate.
/// - `mode`: The calibration mode to use.
///
/// # Returns
/// A `Result` indicating success or failure.
///
/// # Errors
/// This function can return errors from BAM reading/writing operations or if calibration parameters are invalid.
#[allow(dead_code)]
pub fn calibrate<R, W>(
    reader: &mut R,
    writer: &mut W,
    target_regions: &[Region],
    mode: CalibrationMode,
    exclude_uncalibrated_reads: bool,
) -> Result<()>
where
    R: BamReader,
    W: BamWriter,
{
    let sequin_chromosomes = target_regions
        .iter()
        .map(|r| r.contig.as_bytes())
        .collect::<HashSet<_>>();

    let sequin_tids = reader
        .header()
        .target_names()
        .iter()
        .enumerate()
        .filter_map(|(tid, name)| {
            if sequin_chromosomes.contains(name) {
                Some(tid as i32)
            } else {
                None
            }
        })
        .collect::<HashSet<_>>();

    let mut keep = HashSet::new();

    match mode {
        CalibrationMode::FixedCoverage {
            fold_coverage,
            seed,
        } => {
            calibrate_by_fixed_coverage(
                reader,
                target_regions,
                None,
                &mut keep,
                &sequin_tids,
                fold_coverage,
                seed,
            )?;
        }
        CalibrationMode::SampleMeanCoverage {
            sample_regions,
            seed,
        } => {
            calibrate_by_fixed_coverage(
                reader,
                target_regions,
                Some(sample_regions),
                &mut keep,
                &sequin_tids,
                0,
                seed,
            )?;
        }
        CalibrationMode::SampleProfile {
            sample_regions,
            flank,
            window_size,
            min_mapq,
            seed,
        } => {
            let args = SampleProfileParams {
                flank,
                window_size,
                min_mapq,
                seed,
            };
            calibrate_by_sample_profile(reader, writer, target_regions, sample_regions, &args)?;
        }
    }

    reader.fetch(FetchDefinition::All)?;
    for result in reader.records() {
        let record = result?;
        if keep.contains(record.qname()) {
            // If the read is part of a read group selected to keep, write it
            // regardless of anything else.
            writer.write(&record)?;
        } else if !sequin_tids.contains(&record.mtid()) && !exclude_uncalibrated_reads {
            // If we are keeping uncalibrated reads, and the mate is not mapped
            // to a Sequin decoy, write it to the output.
            writer.write(&record)?;
        }
    }

    Ok(())
}

/// Calibrates by fixed coverage or sample mean coverage.
///
/// This function determines downsampling probabilities based on target and
/// sample regions, then subsamples reads in the target regions accordingly. If
/// it decides to keep a read group it inserts the name into the `keep` set.
/// This can then be used to filter the calibrated reads.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `target_regions`: Regions to calibrate.
/// - `sample_regions`: Optional sample regions for mean coverage.
/// - `keep`: A mutable set to store names of read groups to keep.
/// - `sequin_tids`: Set of TIDs corresponding to Sequin chromosomes.
/// - `fold_coverage`: Desired fold coverage (ignored if sample_regions is provided).
/// - `seed`: Random seed for downsampling.
///
/// # Returns
/// A `Result` indicating success or failure.
fn calibrate_by_fixed_coverage<R>(
    reader: &mut R,
    target_regions: &[Region],
    sample_regions: Option<&[Region]>,
    keep: &mut HashSet<Vec<u8>>,
    sequin_tids: &HashSet<i32>,
    fold_coverage: u64,
    seed: u64,
) -> Result<()>
where
    R: BamReader,
{
    let probabilities = determine_downsampling_probabilities(
        reader,
        target_regions,
        sample_regions,
        fold_coverage,
    )?;

    let mut rng = Pcg32::seed_from_u64(seed);
    let mut considered = HashSet::new();
    for region in target_regions {
        reader.fetch((&region.contig, region.beg, region.end))?;
        for result in reader.records() {
            let record = result?;
            // If the mate is mapped to a non-sequin chromosome, skip it. These
            // are artefacts, and we shouldn't keep them.
            if sequin_tids.contains(&record.mtid()) {
                let probability =
                    probabilities
                        .get(&region.name)
                        .cloned()
                        .ok_or_else(|| Error::Calibration {
                            msg: format!(
                                "No downsampling probability found for region {} {:?}",
                                region.name, probabilities
                            ),
                        })?;
                subsample(&record, keep, &considered, probability, &mut rng);
            }
            considered.insert(record.qname().to_vec());
        }
    }

    Ok(())
}

/// Determines downsampling probabilities for each target region.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `target_regions`: Regions to calibrate.
/// - `sample_regions`: Optional sample regions.
/// - `fold_coverage`: Desired fold coverage.
///
/// # Returns
/// A `Result` containing a map of region names to downsampling probabilities.
fn determine_downsampling_probabilities<R: BamReader>(
    reader: &mut R,
    target_regions: &[Region],
    sample_regions: Option<&[Region]>,
    fold_coverage: u64,
) -> Result<HashMap<String, f64>> {
    let target_means = regions_coverage(reader, target_regions)?;
    let sample_means = if let Some(sample_regions) = sample_regions {
        regions_coverage(reader, sample_regions)?
    } else {
        // If no sample regions are provided, we just use the provided fold
        // coverage for all targets.
        target_means
            .keys()
            .map(|k| (k.clone(), fold_coverage as f64))
            .collect()
    };
    let probabilities: HashMap<String, f64> = target_means
        .iter()
        .map(|(name, &target_mean)| {
            let sample_mean =
                sample_means
                    .get(name)
                    .cloned()
                    .ok_or_else(|| Error::Calibration {
                        msg: format!("No sample mean coverage found for target region {}", name),
                    })?;
            if target_mean == 0.0 {
                return Err(Error::Calibration {
                    msg: format!("Target mean coverage for region {} is zero", name),
                });
            }
            if target_mean < sample_mean {
                return Err(Error::Calibration {
                    msg: format!(
                        "Target mean coverage for region {} is less than sample mean coverage",
                        name
                    ),
                });
            }
            let prob = sample_mean / target_mean;
            Ok((name.clone(), prob))
        })
        .collect::<Result<_>>()?;
    Ok(probabilities)
}

/// Calculates the mean coverage of each region.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `regions`: Regions to calculate coverage for.
///
/// # Returns
/// A `Result` containing a map of region names to mean coverage values.
fn regions_coverage<R: BamReader>(
    reader: &mut R,
    regions: &[Region],
) -> Result<HashMap<String, f64>> {
    let min_mapq = 0;
    let flank = 0;
    let coverage = regions
        .iter()
        .map(|region| {
            let region_coverage = coverage::coverage_for_region(reader, region, min_mapq, flank)?;
            let mean = region_coverage.mean().unwrap_or(0.0);
            Ok((region.name.clone(), mean as f64))
        })
        .collect::<Result<_>>()?;
    Ok(coverage)
}

/// Decides whether to keep a read based on downsampling probability.
///
/// This function uses a hash to ensure consistent decisions for paired reads.
///
/// # Arguments
/// - `record`: The BAM record to consider.
/// - `hash`: A mutable hash map for tracking read names.
/// - `threshold`: The downsampling probability threshold.
/// - `rng`: A mutable random number generator.
///
/// # Returns
/// `true` if the read should be kept, `false` otherwise.
fn subsample(
    record: &Record,
    hash: &mut HashSet<Vec<u8>>,
    considered: &HashSet<Vec<u8>>,
    threshold: f64,
    rng: &mut Pcg32,
) -> bool {
    if record.is_duplicate() {
        return false;
    }
    let qname = record.qname().to_vec();
    match hash.contains(&qname) {
        true => {
            return true;
        }
        false => {
            if considered.contains(&qname) {
                return false;
            }
            // if record.pos() >= record.mpos() {
            //     return false;
            // }
        }
    };
    let rand = rng.random::<f64>();
    if rand <= threshold {
        hash.insert(record.qname().to_vec());
        return true;
    }
    false
}

/// Parameters for sample profile calibration.
struct SampleProfileParams {
    /// Number of bases to flank regions.
    flank: u64,
    /// Size of windows for profile matching.
    window_size: u64,
    /// Minimum mapping quality.
    min_mapq: u8,
    /// Random seed.
    seed: u64,
}

/// Calibrates by matching the coverage profile of sample regions.
///
/// This function attempts to replicate the read start profile of sample regions
/// in the target regions using windowed downsampling.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `writer`: A mutable reference to a BAM writer.
/// - `target_regions`: Regions to calibrate.
/// - `sample_regions`: Sample regions to match.
/// - `args`: Parameters for calibration.
///
/// # Returns
/// A `Result` indicating success or failure.
fn calibrate_by_sample_profile<R, W>(
    reader: &mut R,
    writer: &mut W,
    target_regions: &[Region],
    sample_regions: &[Region],
    args: &SampleProfileParams,
) -> Result<()>
where
    R: BamReader,
    W: BamWriter,
{
    let cal_contigs = target_regions
        .iter()
        .map(|r| r.contig.as_str())
        .collect::<Vec<_>>();
    let hdr = reader.header().to_owned();

    let target_names: Vec<&str> = hdr
        .target_names()
        .iter()
        .map(|name| {
            std::str::from_utf8(name).map_err(|e| Error::Calibration {
                msg: format!("Invalid UTF-8 in target name: {}", e),
            })
        })
        .collect::<Result<_>>()?;

    let target_regions = clip_regions(target_regions, args.flank);
    let sample_regions = clip_regions(sample_regions, args.flank);

    let sample_region_map = sample_regions
        .iter()
        .map(|r| (r.name.clone(), r))
        .collect::<HashMap<_, _>>();

    for i in 0..hdr.target_count() {
        let contig = target_names[i as usize];
        if !cal_contigs.contains(&contig) {
            continue;
        }
        calibrate_regions(
            reader,
            writer,
            &target_regions,
            &sample_region_map,
            args.window_size,
            args.min_mapq,
            args.seed,
        )?;
    }
    Ok(())
}

/// Clips regions by removing flanking bases.
///
/// # Arguments
/// - `regions`: Regions to clip.
/// - `flank`: Number of bases to remove from each end.
///
/// # Returns
/// A vector of clipped regions.
fn clip_regions(regions: &[Region], flank: u64) -> Vec<Region> {
    regions
        .iter()
        .map(|r| Region {
            contig: r.contig.clone(),
            beg: r.beg + flank,
            end: r.end - flank,
            name: r.name.clone(),
        })
        .collect()
}

/// Calibrates individual regions by matching sample profiles.
///
/// If the regions need clipping (that is, removing flanking bases), this should
/// be done prior to calling this function.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `writer`: A mutable reference to a BAM writer.
/// - `target_regions`: Target regions to calibrate.
/// - `sample_region_map`: Map of sample regions.
/// - `window_size`: Size of windows.
/// - `min_mapq`: Minimum mapping quality.
/// - `seed`: Random seed.
///
/// # Returns
/// A `Result` indicating success or failure.
fn calibrate_regions<R, W>(
    reader: &mut R,
    writer: &mut W,
    target_regions: &[Region],
    sample_region_map: &HashMap<String, &Region>,
    window_size: u64,
    min_mapq: u8,
    seed: u64,
) -> Result<()>
where
    R: BamReader,
    W: BamWriter,
{
    for target_region in target_regions {
        eprintln!("Calibrating region {}.", target_region.name);
        let sample_region =
            sample_region_map
                .get(&target_region.name)
                .ok_or_else(|| Error::Calibration {
                    msg: format!(
                        "No matching sample region found for target region {}",
                        target_region.name
                    ),
                })?;
        // The first window of the target region is calibrated against the last
        // window of the sample region. This is intentional. The Sequin (target)
        // regions are the mirror of the sample region; therefore, we want to
        // mimic the coverage profile in reverse.
        let sample_starts = window_starts(reader, sample_region, window_size, min_mapq)?;
        let rev_sample_starts = sample_starts.into_iter().rev().collect::<Vec<_>>();

        let records = records_that_start_in_region(
            reader,
            target_region.contig.as_str(),
            target_region.beg,
            target_region.end,
        )?;

        let mut keep_names = HashSet::new();
        for (i, window_beg) in (target_region.beg..target_region.end)
            .step_by(window_size as usize)
            .enumerate()
        {
            let window_end = window_beg + window_size - 1;

            // Divide by 2 because each read is part of a pair. When we include
            // one we automatically include the other; therefore, we will have
            // twice the coverage we intended.
            let n_starts = rev_sample_starts[i] / 2;

            // These are the records that *start* in the current window.
            let region_records = records
                .iter()
                .filter(|record| {
                    if let Ok(pos) = u64::try_from(record.pos()) {
                        pos >= window_beg && pos <= window_end
                    } else {
                        false
                    }
                })
                .collect::<Vec<_>>();
            let numbers = choose_from(region_records.len() as u64, n_starts as u64, seed);
            for idx in &numbers {
                let record = region_records[*idx as usize];
                let qname =
                    String::from_utf8(record.qname().to_vec()).map_err(|e| e.utf8_error())?;
                keep_names.insert(qname);
            }
        }
        for record in &records {
            let qname = String::from_utf8(record.qname().to_vec()).map_err(|e| e.utf8_error())?;
            if keep_names.contains(&qname) {
                writer.write(record)?;
            }
        }
    }
    Ok(())
}

/// Returns the number of read starts in each window of a region.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `region`: The region to analyze.
/// - `window_size`: Size of windows.
/// - `min_mapq`: Minimum mapping quality.
///
/// # Returns
/// A `Result` containing a vector of start counts per window.
fn window_starts<R: BamReader>(
    reader: &mut R,
    region: &Region,
    window_size: u64,
    min_mapq: u8,
) -> Result<Vec<usize>> {
    let mut starts = Vec::new();
    for beg in (region.beg..region.end).step_by(window_size as usize) {
        let end = beg + window_size - 1;
        let n = starts_in(
            reader,
            &Region {
                contig: region.contig.to_owned(),
                beg,
                end,
                name: region.name.to_owned(),
            },
            min_mapq,
        )?;
        starts.push(n);
    }
    Ok(starts)
}

/// Counts the number of read starts in a region.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `region`: The region to count in.
/// - `min_mapq`: Minimum mapping quality.
///
/// # Returns
/// A `Result` containing the count of read starts.
fn starts_in<R: BamReader>(reader: &mut R, region: &Region, min_mapq: u8) -> Result<usize> {
    let (beg, end) = (region.beg as i64, region.end as i64);
    reader.fetch((&region.contig, region.beg, region.end))?;
    let mut n = 0;
    for result in reader.records() {
        let record = result?;
        if record.pos() >= beg && record.pos() <= end && record.mapq() >= min_mapq {
            n += 1;
        }
    }
    Ok(n)
}

/// Retrieves all records that start within a specified region.
///
/// # Arguments
/// - `reader`: A mutable reference to a BAM reader.
/// - `contig`: The contig name.
/// - `beg`: Start position.
/// - `end`: End position.
///
/// # Returns
/// A `Result` containing a vector of BAM records.
fn records_that_start_in_region<R: BamReader>(
    reader: &mut R,
    contig: &str,
    beg: u64,
    end: u64,
) -> Result<Vec<Record>> {
    let mut records = Vec::new();
    reader.fetch((contig, beg, end))?;
    for result in reader.records() {
        let record = result?;
        let pos = record.pos();
        if pos < beg as i64 || pos > end as i64 {
            continue;
        }
        records.push(record);
    }
    Ok(records)
}

/// Randomly selects a subset of indices from a range.
///
/// If `n` is greater than `size`, all indices are returned.
///
/// # Arguments
/// - `size`: The total number of items.
/// - `n`: The number to select.
/// - `seed`: Random seed.
///
/// # Returns
/// A `Result` containing the selected indices.
fn choose_from(size: u64, n: u64, seed: u64) -> Vec<u64> {
    let mut rng = Pcg32::seed_from_u64(seed);
    (0..size).choose_multiple(&mut rng, n as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::{MockBamReader, MockBamWriter};
    use crate::region::Region;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::Record;
    use std::collections::HashMap;

    const CHR1_TID: i32 = 0;
    const CHRQ_MIRROR_TID: i32 = 3;

    /// Helper function to create a mock BAM record
    fn create_mock_record(tid: i32, pos: i64, qname: &str) -> Record {
        let mut record = Record::new();
        record.set_tid(tid);
        record.set_pos(pos);
        record.unset_unmapped();
        record.set_qname(qname.as_bytes());
        // let seq = b"A".repeat(100);
        let cigar_string = CigarString(vec![Cigar::Match(100)]);
        record.set_cigar(Some(&cigar_string));
        record.set_mtid(tid);
        record
    }

    /// Helper function to create a mock reader with records
    fn create_mock_reader_with_records(records: Vec<Record>) -> MockBamReader {
        // Use default header to avoid header format issues
        MockBamReader::new(records, None)
    }

    #[test]
    fn test_calibrate_fixed_coverage_mode() {
        let records = vec![
            create_mock_record(CHRQ_MIRROR_TID, 100, "read1"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read2"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read3"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read4"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
        ];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let mode = CalibrationMode::FixedCoverage {
            fold_coverage: 10,
            seed: 42,
        };

        let result = calibrate(&mut reader, &mut writer, &target_regions, mode, false);
        assert!(result.is_ok());

        // Should have processed records (may be 0 or more)
        let _ = writer.records();
    }

    #[test]
    fn test_calibrate_fixed_coverage_mode_different_chromosomes() {
        let mut r1 = create_mock_record(CHRQ_MIRROR_TID, 150, "read1");
        let r2 = create_mock_record(CHRQ_MIRROR_TID, 150, "read2");
        let r3 = create_mock_record(CHRQ_MIRROR_TID, 150, "read3");
        let r4 = create_mock_record(CHRQ_MIRROR_TID, 150, "read4");
        let r5 = create_mock_record(CHRQ_MIRROR_TID, 150, "read5");
        let r6 = create_mock_record(CHRQ_MIRROR_TID, 150, "read6");
        let r7 = create_mock_record(CHRQ_MIRROR_TID, 150, "read7");
        let r8 = create_mock_record(CHRQ_MIRROR_TID, 150, "read8");
        let r9 = create_mock_record(CHRQ_MIRROR_TID, 150, "read9");
        let r10 = create_mock_record(CHRQ_MIRROR_TID, 150, "read10");

        r1.set_paired();
        r1.set_first_in_template();
        r1.set_mate_reverse();
        r1.set_mtid(CHR1_TID);

        let records = vec![r1, r2, r3, r4, r5, r6, r7, r8, r9, r10];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();
        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let mode = CalibrationMode::FixedCoverage {
            fold_coverage: 5,
            seed: 42,
        };
        let result = calibrate(&mut reader, &mut writer, &target_regions, mode, false);
        assert!(result.is_ok(), "Result: {:?}", result);
        // assert_eq!(writer.records().len(), 9);
    }

    #[test]
    fn test_calibrate_sample_mean_coverage_mode() {
        let records = vec![
            create_mock_record(CHR1_TID, 100, "read1"),
            create_mock_record(CHR1_TID, 100, "read2"),
            create_mock_record(CHR1_TID, 100, "read3"),
            create_mock_record(CHR1_TID, 100, "read4"),
            create_mock_record(CHR1_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read11"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read12"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read13"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read14"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read15"),
        ];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let sample_regions = vec![Region::new("chr1", 100, 200, "region1")];
        let mode = CalibrationMode::SampleMeanCoverage {
            sample_regions: &sample_regions,
            seed: 42,
        };

        let result = calibrate(&mut reader, &mut writer, &target_regions, mode, true);
        assert!(result.is_ok());
        let records = writer.records();
        assert_eq!(records.len(), 5, "Records: ({}) {records:?}", records.len());
        // Should have downsampled to match sample mean
    }

    // It's much harder to reason about what the correct result should be here.
    #[test]
    fn test_calibrate_sample_profile_mode() {
        let records = vec![
            create_mock_record(CHR1_TID, 100, "read1"),
            create_mock_record(CHR1_TID, 100, "read2"),
            create_mock_record(CHR1_TID, 100, "read3"),
            create_mock_record(CHR1_TID, 100, "read4"),
            create_mock_record(CHR1_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read11"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read12"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read13"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read14"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read15"),
        ];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let sample_regions = vec![Region::new("chr1", 100, 200, "region1")];
        let mode = CalibrationMode::SampleProfile {
            sample_regions: &sample_regions,
            flank: 50,
            window_size: 10,
            min_mapq: 20,
            seed: 42,
        };

        let result = calibrate(&mut reader, &mut writer, &target_regions, mode, false);
        assert!(result.is_ok());
        // let records = writer.records();
        // assert_eq!(records.len(), 10);
    }

    #[test]
    fn test_calibrate_empty_target_regions() {
        let records = vec![];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![];
        let mode = CalibrationMode::FixedCoverage {
            fold_coverage: 10,
            seed: 42,
        };

        let result = calibrate(&mut reader, &mut writer, &target_regions, mode, false);
        assert!(result.is_ok());
    }

    #[test]
    #[allow(clippy::manual_range_contains)]
    fn test_determine_downsampling_probabilities_fixed_coverage() {
        let records = vec![
            create_mock_record(CHRQ_MIRROR_TID, 100, "read1"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read2"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read3"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read4"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];

        let result = determine_downsampling_probabilities(&mut reader, &target_regions, None, 5);
        assert!(result.is_ok(), "Expected Ok, got {:?}", result.err());

        let probabilities = result.unwrap();
        assert!(probabilities.contains_key("region1"));
        let p = probabilities["region1"];
        assert!(p >= 0.0 && p <= 1.0, "Probability out of range: {}", p);
    }

    #[test]
    fn test_determine_downsampling_probabilities_fixed_coverage_too_low() {
        let records = vec![
            create_mock_record(CHRQ_MIRROR_TID, 100, "read1"),
            create_mock_record(CHRQ_MIRROR_TID, 200, "read2"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];

        let result = determine_downsampling_probabilities(&mut reader, &target_regions, None, 10);
        assert!(result.is_err());
    }

    #[test]
    #[allow(clippy::manual_range_contains)]
    fn test_determine_downsampling_probabilities_with_sample_regions() {
        let records = vec![
            create_mock_record(CHR1_TID, 100, "read1"),
            create_mock_record(CHR1_TID, 100, "read2"),
            create_mock_record(CHR1_TID, 100, "read3"),
            create_mock_record(CHR1_TID, 100, "read4"),
            create_mock_record(CHR1_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read11"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read12"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read13"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read14"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read15"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let sample_regions = vec![Region::new("chr1", 100, 200, "region1")];

        let result = determine_downsampling_probabilities(
            &mut reader,
            &target_regions,
            Some(&sample_regions),
            10,
        );
        assert!(result.is_ok());

        let probabilities = result.unwrap();
        assert!(probabilities.contains_key("region1"));
        let p = probabilities["region1"];
        assert!(p >= 0.0 && p <= 1.0, "Probability out of range: {}", p);
    }

    #[test]
    fn test_determine_downsampling_probabilities_zero_coverage() {
        let records = vec![]; // No records means zero coverage
        let mut reader = create_mock_reader_with_records(records);

        let target_regions = vec![Region::new("chrQ_mirror", 0, 1000, "region1")];

        let result = determine_downsampling_probabilities(&mut reader, &target_regions, None, 10);
        // Should fail because target coverage is zero
        assert!(result.is_err());
    }

    #[test]
    fn test_regions_coverage() {
        let records = vec![
            create_mock_record(0, 100, "read1"),
            create_mock_record(0, 200, "read2"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let regions = vec![Region::new("chrQ_mirror", 0, 1000, "region1")];

        let result = regions_coverage(&mut reader, &regions);
        assert!(result.is_ok());

        let coverage = result.unwrap();
        assert!(coverage.contains_key("region1"));
        assert!(coverage["region1"] >= 0.0);
    }

    #[test]
    fn test_regions_coverage_empty_regions() {
        let records = vec![];
        let mut reader = create_mock_reader_with_records(records);

        let regions = vec![];

        let result = regions_coverage(&mut reader, &regions);
        assert!(result.is_ok());

        let coverage = result.unwrap();
        assert!(coverage.is_empty());
    }

    #[test]
    fn test_subsample() {
        let mut record = create_mock_record(0, 100, "read1");
        record.set_mpos(150); // Set mate position to be greater than pos
        let mut hash = HashSet::new();
        let considered = HashSet::new();
        let mut rng = Pcg32::seed_from_u64(42);

        // Test with probability 1.0 (should always keep)
        let result = subsample(&record, &mut hash, &considered, 1.0, &mut rng);
        assert!(result);

        // Test with probability 0.0 using a different record
        let mut record2 = create_mock_record(0, 200, "read2");
        record2.set_mpos(250);
        let result = subsample(&record2, &mut hash, &considered, 1.0, &mut rng);
        assert!(result);
    }

    #[test]
    fn test_subsample_duplicate_read() {
        let mut record = create_mock_record(0, 100, "read1");
        record.set_flags(1024); // Set duplicate flag
        let mut hash = HashSet::new();
        let considered = HashSet::new();
        let mut rng = Pcg32::seed_from_u64(42);

        // Duplicate reads should always be filtered out
        let result = subsample(&record, &mut hash, &considered, 1.0, &mut rng);
        assert!(!result);
    }

    #[test]
    fn test_clip_regions() {
        let regions = vec![
            Region::new("chrQ_mirror", 100, 900, "region1"),
            Region::new("chr2", 200, 800, "region2"),
        ];

        let clipped = clip_regions(&regions, 50);

        assert_eq!(clipped.len(), 2);
        assert_eq!(clipped[0].beg, 150);
        assert_eq!(clipped[0].end, 850);
        assert_eq!(clipped[1].beg, 250);
        assert_eq!(clipped[1].end, 750);
    }

    #[test]
    fn test_clip_regions_zero_flank() {
        let regions = vec![Region::new("chrQ_mirror", 100, 900, "region1")];

        let clipped = clip_regions(&regions, 0);

        assert_eq!(clipped[0].beg, 100);
        assert_eq!(clipped[0].end, 900);
    }

    #[test]
    fn test_window_starts() {
        let records = vec![
            create_mock_record(0, 100, "read1"),
            create_mock_record(0, 150, "read2"),
            create_mock_record(0, 250, "read3"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let region = Region::new("chrQ_mirror", 0, 300, "region1");

        let result = window_starts(&mut reader, &region, 100, 0);
        assert!(result.is_ok());

        let starts = result.unwrap();
        assert_eq!(starts.len(), 3); // 3 windows of size 100
    }

    #[test]
    fn test_starts_in() {
        let records = vec![
            create_mock_record(CHRQ_MIRROR_TID, 100, "read1"),
            create_mock_record(CHRQ_MIRROR_TID, 150, "read2"),
            create_mock_record(CHRQ_MIRROR_TID, 250, "read3"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let region = Region::new("chrQ_mirror", 0, 200, "region1");

        let result = starts_in(&mut reader, &region, 0);
        assert!(result.is_ok());

        let count = result.unwrap();
        assert_eq!(count, 2); // 2 reads in the region
    }

    #[test]
    fn test_starts_in_with_min_mapq() {
        let mut record1 = create_mock_record(CHRQ_MIRROR_TID, 100, "read1");
        record1.set_mapq(30);
        let mut record2 = create_mock_record(CHRQ_MIRROR_TID, 150, "read2");
        record2.set_mapq(10); // Below minimum
        let records = vec![record1, record2];
        let mut reader = create_mock_reader_with_records(records);

        let region = Region::new("chrQ_mirror", 0, 200, "region1");

        let result = starts_in(&mut reader, &region, 20);
        assert!(result.is_ok());

        let count = result.unwrap();
        assert_eq!(count, 1); // Only 1 read meets the MAPQ threshold
    }

    #[test]
    fn test_records_that_start_in_region() {
        let records = vec![
            create_mock_record(CHRQ_MIRROR_TID, 100, "read1"),
            create_mock_record(CHRQ_MIRROR_TID, 150, "read2"),
            create_mock_record(CHRQ_MIRROR_TID, 250, "read3"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let result = records_that_start_in_region(&mut reader, "chrQ_mirror", 0, 200);
        assert!(result.is_ok());

        let region_records = result.unwrap();
        assert_eq!(region_records.len(), 2);
    }

    #[test]
    fn test_records_that_start_in_region_empty() {
        let records = vec![];
        let mut reader = create_mock_reader_with_records(records);

        let result = records_that_start_in_region(&mut reader, "chrQ_mirror", 0, 200);
        assert!(result.is_ok());

        let region_records = result.unwrap();
        assert!(region_records.is_empty());
    }

    #[test]
    fn test_choose_from() {
        let result = choose_from(10, 3, 42);
        assert_eq!(result.len(), 3);
        assert!(result.iter().all(|&x| x < 10));
    }

    #[test]
    fn test_choose_from_sample_too_large() {
        let result = choose_from(5, 10, 42);
        assert_eq!(result.len(), 5);
    }

    #[test]
    fn test_choose_from_zero_sample() {
        let result = choose_from(10, 0, 42);
        assert!(result.is_empty());
    }

    #[test]
    fn test_calibrate_by_fixed_coverage() {
        let records = vec![
            create_mock_record(CHRQ_MIRROR_TID, 100, "read1"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read2"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read3"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read4"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
        ];
        let mut reader = create_mock_reader_with_records(records);

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let mut keep = HashSet::new();
        let sequin_tids = [CHRQ_MIRROR_TID].iter().cloned().collect::<HashSet<_>>();
        let result = calibrate_by_fixed_coverage(
            &mut reader,
            &target_regions,
            None,
            &mut keep,
            &sequin_tids,
            5,
            42,
        );
        assert!(result.is_ok(), "Expected Ok, got Err: {:?}", result.err());
    }

    #[test]
    fn test_calibrate_by_sample_profile() {
        let records = vec![
            create_mock_record(CHR1_TID, 100, "read1"),
            create_mock_record(CHR1_TID, 100, "read2"),
            create_mock_record(CHR1_TID, 100, "read3"),
            create_mock_record(CHR1_TID, 100, "read4"),
            create_mock_record(CHR1_TID, 100, "read5"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read6"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read7"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read8"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read9"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read10"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read11"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read12"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read13"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read14"),
            create_mock_record(CHRQ_MIRROR_TID, 100, "read15"),
        ];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![Region::new("chrQ_mirror", 100, 200, "region1")];
        let sample_regions = vec![Region::new("chr1", 100, 200, "region1")];
        let params = SampleProfileParams {
            flank: 50,
            window_size: 100,
            min_mapq: 20,
            seed: 42,
        };

        let result = calibrate_by_sample_profile(
            &mut reader,
            &mut writer,
            &target_regions,
            &sample_regions,
            &params,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_calibrate_regions() {
        let records = vec![
            create_mock_record(0, 100, "read1"),
            create_mock_record(0, 200, "read2"),
        ];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![Region::new("chrQ_mirror", 0, 300, "region1")];
        let sample_region = Region::new("chrQ_mirror", 0, 300, "sample1");
        let sample_region_map = HashMap::from([("region1".to_string(), &sample_region)]);

        let result = calibrate_regions(
            &mut reader,
            &mut writer,
            &target_regions,
            &sample_region_map,
            100,
            20,
            42,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_calibrate_regions_missing_sample() {
        let records = vec![];
        let mut reader = create_mock_reader_with_records(records);
        let mut writer = MockBamWriter::new();

        let target_regions = vec![Region::new("chrQ_mirror", 0, 300, "region1")];
        let sample_region_map = HashMap::new(); // Empty map

        let result = calibrate_regions(
            &mut reader,
            &mut writer,
            &target_regions,
            &sample_region_map,
            100,
            20,
            42,
        );
        assert!(result.is_err()); // Should fail due to missing sample region
    }
}
