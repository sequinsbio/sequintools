/// Provides a trait and implementation for reading BAM files using the
/// `rust-htslib` crate, enabling access to BAM records, headers, and
/// region-based queries with multi-threading support.
///
/// # Overview
/// - Defines the `BamReader` trait, which abstracts over BAM file reading
///   functionality.
/// - Implements `BamReader` for `HtslibBamReader`, a wrapper around
///   `rust_htslib::bam::IndexedReader`.
/// - Includes a mock implementation for testing purposes.
///
/// # Features
/// - Iteration over BAM records via an associated iterator type.
/// - Fetching records from specific regions using `fetch`.
/// - Access to BAM header information.
/// - Multi-threaded reading via `set_threads`.
/// - Reference genome specification with `set_reference`.
///
/// # Usage
/// Implement the `BamReader` trait for custom BAM readers, or use the provided
/// `HtslibBamReader` for standard indexed BAM file access.
///
/// # Testing
/// Includes a `MockBamReader` for unit testing code that depends on the
/// `BamReader` trait.
///
/// # Dependencies
/// - [rust-htslib](https://docs.rs/rust-htslib)
///
/// # Example
/// ```ignore
/// use crate::bam::{BamReader, HtslibBamReader};
/// let reader = HtslibBamReader::from_path(&bam_path)?;
/// for record in reader.records() {
///     // process record
/// }
/// ```
use crate::errors::Result;
use rust_htslib::bam::{FetchDefinition, HeaderView, IndexedReader, Read, Record};
use std::path::Path;

/// A trait for reading BAM files, providing an interface for accessing records,
/// headers, and controlling reading behavior.
///
/// Types implementing this trait can iterate over records, fetch specific
/// regions, set the number of threads for reading, and specify a reference
/// genome.
///
/// # Associated Types
/// - `RecordsIter<'a>`: An iterator over `Result<Record, rust_htslib::errors::Error>`.
///
/// # Required Methods
/// - `header(&self) -> &HeaderView`: Returns a reference to the BAM header.
/// - `fetch(&mut self, definition)`: Restricts the reader to a specific region or definition.
/// - `records(&mut self)`: Returns an iterator over records in the current region or file.
/// - `set_threads(&mut self, n)`: Sets the number of threads for reading.
/// - `set_reference(&mut self, reference)`: Sets the reference genome for the reader.
pub trait BamReader {
    /// Associated iterator type for records
    type RecordsIter<'a>: Iterator<Item = std::result::Result<Record, rust_htslib::errors::Error>>
    where
        Self: 'a;

    /// Returns a reference to the BAM header
    fn header(&self) -> &HeaderView;

    /// Fetches a region or reference specified by the definition
    fn fetch<'a, T: Into<FetchDefinition<'a>>>(
        &mut self,
        definition: T,
    ) -> std::result::Result<(), rust_htslib::errors::Error>;

    /// Returns an iterator over BAM records
    fn records(&mut self) -> Self::RecordsIter<'_>;

    /// Sets the number of threads to use for reading
    fn set_threads(&mut self, n: usize) -> std::result::Result<(), rust_htslib::errors::Error>;

    /// Sets the reference FASTA file for CRAM/BAM decoding
    fn set_reference<P: AsRef<Path>>(
        &mut self,
        reference: P,
    ) -> std::result::Result<(), rust_htslib::errors::Error>;
}

/// BAM/CRAM reader implementation using rust-htslib's IndexedReader.
pub struct HtslibBamReader {
    reader: IndexedReader,
}

impl HtslibBamReader {
    /// Create a new HtslibBamReader from a file path
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = IndexedReader::from_path(path)?;
        Ok(Self { reader })
    }
}

impl BamReader for HtslibBamReader {
    type RecordsIter<'a> = rust_htslib::bam::Records<'a, IndexedReader>;

    /// Returns a reference to the BAM header
    fn header(&self) -> &HeaderView {
        self.reader.header()
    }

    /// Fetches a region or reference specified by the definition
    fn fetch<'a, T: Into<FetchDefinition<'a>>>(
        &mut self,
        definition: T,
    ) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.reader.fetch(definition)
    }

    /// Returns an iterator over BAM records
    fn records(&mut self) -> Self::RecordsIter<'_> {
        self.reader.records()
    }

    /// Sets the number of threads to use for reading
    fn set_threads(&mut self, n: usize) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.reader.set_threads(n)
    }

    /// Sets the reference FASTA file for CRAM/BAM decoding
    fn set_reference<P: AsRef<Path>>(
        &mut self,
        reference: P,
    ) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.reader.set_reference(reference)
    }
}

/// Mock implementation of `BamReader` for testing purposes.
///
/// This struct simulates the behavior of a BAM/CRAM reader and is intended for
/// use in unit tests where reading from actual files is not desired. It allows
/// you to inject a vector of `Record` objects and provides the same trait
/// interface as a real BAM reader.
///
/// # Example
/// ```ignore
/// use crate::bam::{MockBamReader, BamReader};
/// use rust_htslib::bam::Record;
///
/// // Create some mock records
/// let records = vec![Record::new(), Record::new()];
/// let mut reader = MockBamReader::new(records);
///
/// // Iterate over records as you would with a real BAM reader
/// for result in reader.records() {
///     let record = result.unwrap();
///     // ... test logic ...
/// }
/// ```
#[cfg(test)]
pub struct MockBamReader {
    records: Vec<Record>,
    header: HeaderView,
}

#[cfg(test)]
impl MockBamReader {
    /// Create a new MockBamReader from a vector of records
    pub fn new(records: Vec<Record>) -> Self {
        let mut header = rust_htslib::bam::Header::new();

        header.push_record(&rust_htslib::bam::header::HeaderRecord::new(
            b"SQ\tSN:chrQ_mirror\tLN:83800\n",
        ));
        let header_view = HeaderView::from_header(&header);
        Self {
            records,
            header: header_view,
        }
    }
}

#[cfg(test)]
pub struct MockRecords<'a> {
    records: std::slice::Iter<'a, Record>,
}

#[cfg(test)]
impl<'a> Iterator for MockRecords<'a> {
    type Item = std::result::Result<Record, rust_htslib::errors::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next().map(|r| Ok(r.clone()))
    }
}

#[cfg(test)]
impl BamReader for MockBamReader {
    type RecordsIter<'a> = MockRecords<'a>;

    /// Returns a reference to the BAM header (not implemented for mock)
    fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Mock fetch does nothing and always returns Ok
    fn fetch<'a, T: Into<FetchDefinition<'a>>>(
        &mut self,
        _definition: T,
    ) -> std::result::Result<(), rust_htslib::errors::Error> {
        Ok(())
    }

    /// Returns an iterator over the mock records
    fn records(&mut self) -> Self::RecordsIter<'_> {
        MockRecords {
            records: self.records.iter(),
        }
    }

    /// Mock set_threads does nothing and always returns Ok
    fn set_threads(&mut self, _n: usize) -> std::result::Result<(), rust_htslib::errors::Error> {
        Ok(())
    }

    /// Mock set_reference does nothing and always returns Ok
    fn set_reference<P: AsRef<Path>>(
        &mut self,
        _reference: P,
    ) -> std::result::Result<(), rust_htslib::errors::Error> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_htslib_bam_reader_creation() {
        // This test requires an actual BAM file, so we'll test the error case
        let result = HtslibBamReader::from_path(PathBuf::from("nonexistent.bam"));
        assert!(result.is_err());
    }

    #[test]
    fn test_mock_bam_reader_basic() {
        let records = vec![];
        let mut mock_reader = MockBamReader::new(records);

        // Test that records iterator works (should be empty)
        let collected: Vec<_> = mock_reader.records().collect();
        assert_eq!(collected.len(), 0);
    }

    #[test]
    fn test_mock_bam_reader_with_records() {
        // Create a simple record for testing
        let mut record = Record::new();
        record.set_tid(0);
        record.set_pos(100);

        let records = vec![record];
        let mut mock_reader = MockBamReader::new(records);

        // Test that we get back our record
        let collected: Vec<_> = mock_reader.records().collect();
        assert_eq!(collected.len(), 1);
        assert!(collected[0].is_ok());
    }

    /// Example of testing a function that takes a BamReader
    #[test]
    fn test_example_function_with_mock() {
        // This shows how to test functions that accept BamReader trait objects
        fn process_bam_records<B: BamReader>(mut bam_reader: B) -> usize {
            bam_reader.records().count()
        }

        let records = vec![Record::new(), Record::new(), Record::new()];
        let mock_reader = MockBamReader::new(records);

        let count = process_bam_records(mock_reader);
        assert_eq!(count, 3);
    }

    #[test]
    fn test_mock_bam_reader_methods_dont_error() {
        let mut mock_reader = MockBamReader::new(vec![]);

        // All these methods should succeed on mock
        assert!(mock_reader.set_threads(4).is_ok());
        assert!(mock_reader.set_reference("/path/to/ref.fa").is_ok());
        assert!(mock_reader.fetch((0, 100, 200)).is_ok());
    }

    // This exists because codecov doesn't like the fact that the real
    // implementation isn't tested even though that's the point of the mock.
    #[test]
    fn test_htslib_bam_reader() {
        let bam_path = PathBuf::from("testdata/calibrated.bam");
        let mut reader = HtslibBamReader::from_path(&bam_path).expect("Failed to open BAM file");

        assert!(reader.set_threads(2).is_ok(), "Failed to set threads");

        let header = reader.header();
        assert!(header.target_count() > 0, "Header has no targets");

        let reference = PathBuf::from("testdata/reference.fasta");
        assert!(
            reader.set_reference(reference).is_ok(),
            "Failed to set reference"
        );

        reader
            .fetch(FetchDefinition::All)
            .expect("Failed to fetch all");
        let n = reader.records().count();
        let expected = 13321;
        assert_eq!(
            n, expected,
            "Unexpected number of records: expected {expected}, got {n}"
        );
    }
}
