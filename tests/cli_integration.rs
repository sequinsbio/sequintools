use std::fs;
use std::fs::File;
use std::process::{Command, Stdio};
use tempfile::TempDir;

#[test]
fn test_cli_bedcov() {
    let expected_output = "\
name,chrom,beg,end,min,max,mean,std,cv
variant_1,chrQ_mirror,200,3200,0,51,34.56,11.56,0.33
variant_2,chrQ_mirror,3400,6400,0,50,30.89,11.59,0.38
variant_3,chrQ_mirror,6600,9600,0,57,36.69,12.90,0.35";
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
    let expected_md5 = "a7c2c35e7fe0f9459c96d93cc893aa70";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
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
    let expected_md5 = "524db272c7eb68049a3c4b72ab307f35";
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
    let expected_md5 = "8071d48d8130b2b6a7d459b51b662356";
    assert_eq!(computed_md5, expected_md5, "MD5 checksum does not match");
}
