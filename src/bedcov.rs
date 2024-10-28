use crate::calibrate::mean_depth;
use crate::read_bed;
use anyhow::Result;
use rust_htslib::bam;
use std::fs::File;
use std::io;

// For each region in `bed_path`, calculate the mean depth. This should be using
// the mean_depth function in calibrate.
pub fn bedcov(bam_path: String, bed_path: String, min_mapq: u8, flank: i32) -> Result<()> {
    let f = File::open(bed_path)?;
    let mut reader = io::BufReader::new(f);
    let regions = read_bed(&mut reader)?;
    let mut bam = bam::IndexedReader::from_path(bam_path).unwrap();
    let mut wtr = csv::Writer::from_writer(io::stdout());
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
