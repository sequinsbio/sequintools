use rust_htslib::bam::{self, Read};
use std::fs;
use std::fs::File;
use std::path::Path;
use std::process::{Command, Stdio};
use tempfile::TempDir;

fn calculate_md5_without_pg_records(bam_path: &Path) -> String {
    let mut samtools_output = Command::new("samtools")
        .args(["view", "-h", bam_path.to_str().unwrap()])
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to execute samtools command");
    let grep_output = Command::new("grep")
        .args(["-v", "^@PG"])
        .stdin(samtools_output.stdout.take().unwrap())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to execute grep command");
    let output = &grep_output
        .wait_with_output()
        .expect("Failed to read grep output");
    samtools_output.wait().expect("Samtools command failed");
    if !output.status.success() {
        panic!(
            "grep command failed with status: {}\nstderr: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let digest = md5::compute(&output.stdout);
    format!("{digest:x}")
}

#[test]
fn test_cli_bedcov() {
    let expected_output = "\
name,chrom,beg,end,min,max,mean,std,cv
variant_1,chrQ_mirror,200,3200,0,172,116.84,39.77,0.34
variant_2,chrQ_mirror,3400,6400,0,163,113.95,37.87,0.33
variant_3,chrQ_mirror,6600,9600,0,163,116.94,39.02,0.33";

    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.csv");
    let file = File::create(&output_path).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "bedcov",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "testdata/calibrated.bam",
        ])
        .stdout(Stdio::from(file))
        .stderr(Stdio::piped())
        .output()
        .expect("Failed to execute command");
    assert!(
        output.status.success(),
        "Command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output_path.exists());
    let contents = fs::read_to_string(output_path).unwrap();
    assert_eq!(contents.trim(), expected_output.trim());
}

#[test]
fn test_calibrate_fixed_coverage() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "calibrate",
            "--bed",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "-o",
            output_path.to_str().unwrap(),
            "--write-index",
            "testdata/uncalibrated.bam",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let computed_md5 = calculate_md5_without_pg_records(&output_path);
    let expected_md5 = "3c9bbf17b331d6bf1ddcc1f5725175ec";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
    let mut index_path = output_path.clone();
    index_path.set_extension("bam.bai");
    assert!(index_path.exists(), "No such file: {index_path:?}");
}

#[test]
fn test_calibrate_sample_mean_coverage() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "calibrate",
            "--bed",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "--sample-bed",
            "testdata/resources/sequin_regions.hg38.bed",
            "-o",
            output_path.to_str().unwrap(),
            "testdata/uncalibrated.bam",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let computed_md5 = calculate_md5_without_pg_records(&output_path);
    let expected_md5 = "967f061c503316128bcccd605c435224";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
}

#[test]
fn test_calibrate_sample_coverage_profile() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "calibrate",
            "--experimental",
            "--bed",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "--sample-bed",
            "testdata/resources/sequin_regions.hg38.bed",
            "-o",
            output_path.to_str().unwrap(),
            "testdata/uncalibrated.bam",
        ])
        .output()
        .expect("Failed to execute command");

    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let computed_md5 = calculate_md5_without_pg_records(&output_path);
    let expected_md5 = "d399ae99545464c9631ace3258e57860";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
}

#[test]
fn test_calibrate_cram_input() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "calibrate",
            "--bed",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "-o",
            output_path.to_str().unwrap(),
            "-T",
            "testdata/genome_with_sequins.fasta",
            "testdata/uncalibrated.cram",
        ])
        .output()
        .expect("Failed to execute command");
    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let computed_md5 = calculate_md5_without_pg_records(&output_path);
    let expected_md5 = "bae7f889aa8a1c1f7db5d10d2b2dd0ff";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
}

// We can't check MD5 because the CRAM output may vary in non-deterministic ways
// on different machines. This is due to the choice of compression level and
// other factors. However, all output files (BAM and CRAM) should have the same
// content.

// These aren't very robust tests, but it's difficult to check every aspect of
// every read to ensure identity.

#[test]
fn test_calibrate_cram_output() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.cram");
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "calibrate",
            "--bed",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "-T",
            "testdata/genome_with_sequins.fasta",
            "--cram",
            "-o",
            output_path.to_str().unwrap(),
            "--write-index",
            "testdata/uncalibrated.bam",
        ])
        .output()
        .expect("Failed to execute command");
    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut reader =
        bam::Reader::from_path(&output_path).expect("Should be able to open output CRAM");
    let record_count = reader.records().count();
    assert_eq!(record_count, 4718);
    let mut index_path = output_path.clone();
    index_path.set_extension("cram.crai");
    assert!(index_path.exists(), "No such file: {index_path:?}");
}

#[test]
fn test_calibrate_cram_input_output() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.cram");
    let output = Command::new(env!("CARGO_BIN_EXE_sequintools"))
        .args([
            "calibrate",
            "--bed",
            "testdata/resources/sequin_regions.chrQ_mirror.bed",
            "-T",
            "testdata/genome_with_sequins.fasta",
            "--cram",
            "-o",
            output_path.to_str().unwrap(),
            "testdata/uncalibrated.cram",
        ])
        .output()
        .expect("Failed to execute command");
    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut reader =
        bam::Reader::from_path(output_path).expect("Should be able to open output CRAM");
    let record_count = reader.records().count();
    assert_eq!(record_count, 4718);
}
