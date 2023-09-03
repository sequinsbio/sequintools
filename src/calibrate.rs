use crate::read_bed;
use crate::CalibrateArgs;
use crate::Region;
use rand::seq::IteratorRandom;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::fmt;
use std::process::exit;
use std::vec::Vec;

pub fn calibrate(args: CalibrateArgs) -> Result<(), Box<dyn Error>> {
    let _ = match args.sample_bed {
        Some(_) => calibrate_by_sample_coverage(args),
        None => calibrate_by_standard_coverage(args),
    };
    Ok(())
}

// fn calibrate_by_mean_sample_coverage() {}

// pub fn calibrate(flank: u64, seed: u64, output: String, bed: String, path: String) {
pub fn calibrate_by_standard_coverage(args: CalibrateArgs) -> Result<(), Box<dyn Error>> {
    let regions = read_bed(&args.bed);

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
            let depth = mean_depth(&mut bam, &region, args.flank);
            println!(
                "{}:{}-{} {} {}",
                region.contig,
                region.beg,
                region.end,
                depth,
                depth * (coverage / depth)
            );

            let chrom = region.contig.as_str();
            let beg = region.beg + args.flank;
            let end = region.end - args.flank;

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
        match &args.output {
            Some(path) => bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap(),
            None => (), //
        }
    }
    Ok(())
}

fn mean_depth(bam: &mut bam::IndexedReader, region: &Region, flank: u64) -> f64 {
    let chrom = region.contig.as_str();
    let beg = region.beg + flank;
    let end = region.end - flank;
    bam.fetch((chrom, beg, end)).unwrap();
    let mut n: u64 = 0;
    let mut c: u64 = 0;
    let include_flags: u16 = 0;
    let exclude_flags: u16 = 0;
    let min_mapq: u8 = 10;
    for p in bam.pileup() {
        let pileup = p.unwrap();
        let pileup_pos = pileup.pos() as u64;
        if pileup_pos < beg || pileup_pos > end {
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
            .count();

        n += 1;
        c += u64::from(depth as u64);
    }
    c as f64 / n as f64
}

//  What if we use a sliding window over the region and determine the depth of
//  that window in the sample data. For the sequin data we randomly sample but
//  only if the read starts in that window. We keep both ends of a pair. How do
//  we handle secondary/supplementry/duplicate reads? We can experiment with
//  different window sizes to get the best replication of the sample coverage.
fn calibrate_by_sample_coverage(args: CalibrateArgs) -> Result<(), Box<dyn Error>> {
    // shouldn't be calling this function if sample_bed is not provided
    let sample_bed = args.sample_bed.as_ref().unwrap();
    let cal_contigs = calibrated_contigs(&args.bed);
    let regions = read_bed(&args.bed);
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
            Some(path) => bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap(),
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
            if cal_contigs.contains(&contig.to_string()) {
                let mut sample_regions_map = HashMap::new();
                let sample_regions = read_bed(&sample_bed);
                for region in &sample_regions {
                    sample_regions_map.insert(&region.name, region);
                }
                calibrate_regions(&mut bam, &mut out, &regions, sample_regions_map, &args);
            } else {
                if !args.exclude_uncalibrated_reads {
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
    }
    if args.write_index {
        match &args.output {
            Some(path) => {
                eprintln!("Writing index");
                bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
            }
            None => (), // can't index stdout
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
            args.flank,
            args.window_size,
            args.min_mapq,
        );
        let rev_sample_starts: Vec<i64> = sample_starts.into_iter().rev().collect();

        let region_beg = region.beg + args.flank;
        let region_end = region.end - args.flank - args.window_size;
        // let region_beg = region.beg;
        // let region_end = region.end - args.window_size;

        // This should be all reads mapped to a region
        let records = records_that_start_in_region(bam, contig, region.beg, region.end);
        let mut keep_names = HashSet::new();
        let mut i = 0;

        for window_beg in (region_beg..region_end).step_by(args.window_size as usize) {
            let window_end = window_beg + args.window_size - 1;

            // We are always keeping both reads from paired-end data, so
            // everytime we select a read to keep we are actually keeping two
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
            i += 1;
        }

        for record in &records {
            let qname = String::from_utf8(record.qname().to_vec()).unwrap(); // really? there must be a better way
            if keep_names.contains(&qname) {
                out.write(&record).unwrap();
            }
        }
    }
}

fn starts_in(bam: &mut bam::IndexedReader, region: &Region, min_mapq: u8) -> i64 {
    // let min_mapq = 10;
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
    let region_end = region.end - flank - window_size as u64;
    // let region_beg = region.beg;
    // let region_end = region.end - window_size as u64;
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
    // let mut xs = Vec::new();
    // for _ in 0..n {
    //     let i = (0..size).choose(&mut rng).unwrap();
    //     xs.push(i);
    // }
    // Ok(xs)
}

// Return the unique set of contig names in a BED file
fn calibrated_contigs(path: &String) -> HashSet<String> {
    let mut contigs = HashSet::new();
    let regions = read_bed(path);
    for region in &regions {
        let contig = String::from(&region.contig);
        contigs.insert(contig);
    }
    contigs
}
