use crate::read_bed;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io;

// For each region in `bed_path`, calculate the mean depth
pub fn bedcov(bam_path: String, bed_path: String, min_mapq: u8) -> Result<()> {
    let f = File::open(&bed_path)?;
    let mut reader = io::BufReader::new(f);
    let regions = read_bed(&mut reader)?;
    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();

    let include_flags: u16 = 0;
    let exclude_flags: u16 = 0;

    for region in &regions {
        let mut n = 0;
        // let mut c = 0;
        let mut o = 0;
        let chrom = region.contig.as_str();
        let beg = region.beg;
        let end = region.end;
        bam.fetch((chrom, beg, end)).unwrap();
        for p in bam.pileup() {
            let pileup = p.unwrap();

            let pileup_pos = pileup.pos() as u64;
            if pileup_pos < beg || pileup_pos > end {
                continue;
            }

            // let depth = pileup.depth();

            let other_depth = pileup
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
            // c += depth as i64;
            o += other_depth as i64;
        }
        // let mean_depth = c as f64 / n as f64;
        let other_mean_depth = o as f64 / n as f64;
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.2}",
            chrom, beg, end, region.name, o, other_mean_depth
        );
    }
    Ok(())
}
