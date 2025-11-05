use rust_htslib::bam::{self, Read};
use std::fs;
use std::fs::File;
use std::process::{Command, Stdio};
use tempfile::TempDir;

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
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
            "--",
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
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
            "--",
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
    let file_bytes = fs::read(&output_path).expect("Should have been able to read the file");
    let digest = md5::compute(&file_bytes);
    let computed_md5 = format!("{digest:x}");
    let expected_md5 = "f7441d6fac943fcc1649af48d34ef77c";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
    let mut index_path = output_path.clone();
    index_path.set_extension("bam.bai");
    assert!(index_path.exists(), "No such file: {index_path:?}");
}

#[test]
fn test_calibrate_sample_mean_coverage() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
            "--",
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
    let file_bytes = fs::read(&output_path).expect("Should have been able to read the file");
    let digest = md5::compute(&file_bytes);
    let computed_md5 = format!("{digest:x}");
    let expected_md5 = "c074c15845ef58ef9ade8498e88ef264";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
}

#[test]
fn test_calibrate_sample_coverage_profile() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
            "--",
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
    let file_bytes = fs::read(&output_path).expect("Should have been able to read the file");
    let digest = md5::compute(&file_bytes);
    let computed_md5 = format!("{digest:x}");
    let expected_md5 = "92b7f0546a12e187afb5795d4e654f78";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
}

#[test]
fn test_calibrate_cram_input() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("calibrated.bam");
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
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
    let file_bytes = fs::read(&output_path).expect("Should have been able to read the file");
    let digest = md5::compute(&file_bytes);
    let computed_md5 = format!("{digest:x}");
    let expected_md5 = "64bf98716f2bc26ab82e3b4a0be24b55";
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
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
            "--",
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
    let output = Command::new("cargo")
        .args([
            "run",
            "--bin",
            "sequintools",
            "--",
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
