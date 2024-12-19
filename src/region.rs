//! # Region Module
//!
//! This module provides functionality for representing and handling genomic regions within contigs.
//!
//! ## Structs
//!
//! - `Region`: Represents a genomic region with a contig name, start and end positions, and an identifier.
//!
//! ## Implementations
//!
//! - Implements the `Display` trait for `Region` to enable formatted string representations.
//!
//! ## Functions
//!
//! - `load_from_bed`: Loads genomic regions from a BED file, parsing each line into a `Region` struct.
//!
//! ## Tests
//!
//! Contains unit tests for verifying the functionality of the `Region` struct and the `load_from_bed` function.
use anyhow::{bail, Context, Result};
use std::fmt;
use std::io::Read;

/// Represents a genomic region within a contig.
///
/// # Fields
/// - `contig`: The contig where the region is located.
/// - `beg`: The starting position of the region.
/// - `end`: The ending position of the region.
/// - `name`: The identifier name of the region.
#[derive(Debug, PartialEq, Eq)]
pub struct Region {
    pub contig: String,
    pub beg: u64,
    pub end: u64,
    pub name: String,
}

/// Impl Display for Region
impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.beg, self.end)
    }
}

/// Loads genomic regions from a BED file.
///
/// Parses each line of the BED file and constructs a `Region` struct for each valid entry.
///
/// # Arguments
///
/// * `reader` - A mutable reference to a type that implements the `Read` trait, typically a file or buffer.
///
/// # Returns
///
/// A `Result` containing a vector of `Region` structs if successful, or an error if parsing fails.
pub fn load_from_bed<R: Read>(reader: &mut R) -> Result<Vec<Region>> {
    let mut result = Vec::new();
    let mut contents = String::new();
    reader.read_to_string(&mut contents)?;
    for (i, line) in contents.lines().enumerate() {
        let bits: Vec<&str> = line.split_whitespace().collect();
        let [contig, beg_str, end_str, name, ..] = bits[..] else {
            bail!(
                "Incorrect number of columns detected, expected >= 4 found {} (line = {})",
                bits.len(),
                i + 1
            )
        };

        let beg: u64 = beg_str.parse().with_context(|| {
            format!(
                "Beg column is not an integer: is {} (line = {})",
                bits[1],
                i + 1
            )
        })?;
        let end: u64 = end_str.parse().with_context(|| {
            format!(
                "End column is not an integer: is {} (line = {})",
                bits[2],
                i + 1
            )
        })?;
        result.push(Region {
            contig: contig.to_owned(),
            beg,
            end,
            name: name.to_owned(),
        });
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{self, Cursor};

    #[test]
    fn display_region() {
        let region = Region {
            contig: "chr1".to_owned(),
            beg: 100,
            end: 200,
            name: "test_region".to_owned(),
        };
        assert_eq!(region.to_string(), "chr1:100-200");

        let region2 = Region {
            contig: "chrX".to_owned(),
            beg: 0,
            end: 1000,
            name: "test_region2".to_owned(),
        };
        assert_eq!(region2.to_string(), "chrX:0-1000");
    }

    #[test]
    fn load_bed_to_vec() {
        let data = b"chr1\t1\t10\treg1";
        let mut cursor = Cursor::new(data);
        let result = load_from_bed(&mut cursor).unwrap();
        assert_eq!(
            result,
            vec![Region {
                contig: "chr1".to_owned(),
                beg: 1,
                end: 10,
                name: "reg1".to_owned()
            }]
        );
    }

    #[test]
    fn load_multi_lines_bed() {
        let mut cursor = Cursor::new(b"chr1\t1\t10\treg1\nchr2\t2\t20\treg2");
        let result = load_from_bed(&mut cursor).unwrap();
        assert_eq!(
            result,
            vec![
                Region {
                    contig: "chr1".to_owned(),
                    beg: 1,
                    end: 10,
                    name: "reg1".to_owned()
                },
                Region {
                    contig: "chr2".to_owned(),
                    beg: 2,
                    end: 20,
                    name: "reg2".to_owned()
                }
            ]
        );
    }

    #[test]
    fn load_with_bad_start() {
        let mut cursor = Cursor::new(b"chr1\txxx\t10\treg1");
        let err = load_from_bed(&mut cursor).unwrap_err();
        // check that the error message is correct
        assert!(err
            .to_string()
            .contains("Beg column is not an integer: is xxx (line = 1)"));
    }

    #[test]
    fn load_with_bad_end() {
        let mut cursor = Cursor::new(b"chr1\t1\txxx\treg1");
        let err = load_from_bed(&mut cursor).unwrap_err();
        // check that the error message is correct
        assert!(err
            .to_string()
            .contains("End column is not an integer: is xxx (line = 1)"));
    }

    #[test]
    fn load_without_name() {
        let mut cursor = Cursor::new(b"chr1\t1\t10");
        let err = load_from_bed(&mut cursor).unwrap_err();
        // check that the error message is correct
        assert!(err
            .to_string()
            .contains("Incorrect number of columns detected, expected >= 4 found 3"));
    }

    struct ErrorReader;
    impl Read for ErrorReader {
        fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::Other, "bad read"))
        }
    }
    #[test]
    fn load_invalid_reader() {
        let mut reader = ErrorReader;
        let result = load_from_bed(&mut reader);
        assert!(result.is_err());
    }
}
