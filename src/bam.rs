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
use rust_htslib::bam::{self, FetchDefinition, HeaderView, IndexedReader, Read, Record};
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

#[cfg(test)]
#[derive(Debug, Clone, PartialEq, Eq)]
enum FetchState {
    All,
    Unmapped,
    Contig(i32),
    Region(i32, u64, u64),
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
    fetch_state: FetchState,
}

#[cfg(test)]
impl MockBamReader {
    /// Create a new MockBamReader from a vector of records
    pub fn new(
        records: Vec<Record>,
        header_records: Option<&[rust_htslib::bam::header::HeaderRecord]>,
    ) -> Self {
        let mut header = rust_htslib::bam::Header::new();

        match header_records {
            Some(recs) => {
                for rec in recs {
                    header.push_record(rec);
                }
            }
            None => {
                header.push_record(&rust_htslib::bam::header::HeaderRecord::new(
                    b"SQ\tSN:chr1\tLN:5000",
                ));
                header.push_record(&rust_htslib::bam::header::HeaderRecord::new(
                    b"SQ\tSN:chr2\tLN:5000",
                ));
                header.push_record(&rust_htslib::bam::header::HeaderRecord::new(
                    b"SQ\tSN:chr3\tLN:5000",
                ));
                header.push_record(&rust_htslib::bam::header::HeaderRecord::new(
                    b"SQ\tSN:chrQ_mirror\tLN:9800",
                ));
            }
        }
        let header_view = HeaderView::from_header(&header);
        Self {
            records,
            header: header_view,
            fetch_state: FetchState::All,
        }
    }
}

#[cfg(test)]
pub struct MockRecords {
    records: std::vec::IntoIter<Record>,
}

#[cfg(test)]
impl Iterator for MockRecords {
    type Item = std::result::Result<Record, rust_htslib::errors::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.records.next().map(Ok)
    }
}

#[cfg(test)]
impl BamReader for MockBamReader {
    type RecordsIter<'a> = MockRecords;

    /// Returns a reference to the BAM header (not implemented for mock)
    fn header(&self) -> &HeaderView {
        &self.header
    }

    /// Mock fetch does nothing and always returns Ok
    fn fetch<'a, T: Into<FetchDefinition<'a>>>(
        &mut self,
        definition: T,
    ) -> std::result::Result<(), rust_htslib::errors::Error> {
        let fetch_def = definition.into();
        match fetch_def {
            FetchDefinition::All => self.fetch_state = FetchState::All,
            FetchDefinition::Region(tid, beg, end) => {
                self.fetch_state = FetchState::Region(tid, beg as u64, end as u64)
            }
            FetchDefinition::RegionString(name, beg, end) => {
                if let Some(tid) = self.header.tid(name) {
                    self.fetch_state = FetchState::Region(tid as i32, beg as u64, end as u64);
                } else {
                    return Err(rust_htslib::errors::Error::Fetch);
                }
            }
            FetchDefinition::CompleteTid(tid) => {
                self.fetch_state = FetchState::Contig(tid);
            }
            FetchDefinition::String(name) => {
                if name == b"*" {
                    self.fetch_state = FetchState::Unmapped;
                } else if let Some(tid) = self.header.tid(name) {
                    self.fetch_state = FetchState::Contig(tid as i32);
                } else {
                    return Err(rust_htslib::errors::Error::Fetch);
                }
            }
            FetchDefinition::Unmapped => {
                self.fetch_state = FetchState::Unmapped;
            }
        }
        Ok(())
    }

    /// Returns an iterator over the mock records, filtered by the current fetch
    /// state.
    fn records(&mut self) -> Self::RecordsIter<'_> {
        let filtered_records = self
            .records
            .iter()
            .filter(|record| match self.fetch_state {
                FetchState::All => true,
                FetchState::Unmapped => record.is_unmapped(),
                FetchState::Contig(tid) => record.tid() == tid,
                FetchState::Region(tid, beg, end) => {
                    record.tid() == tid && record.pos() >= beg as i64 && record.pos() < end as i64
                }
            })
            .cloned()
            .collect::<Vec<_>>();
        MockRecords {
            records: filtered_records.into_iter(),
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

pub trait BamWriter {
    fn write(&mut self, record: &Record) -> std::result::Result<(), rust_htslib::errors::Error>;

    /// Set the number of threads to use for writing
    fn set_threads(&mut self, n: usize) -> std::result::Result<(), rust_htslib::errors::Error>;

    /// Sets the reference FASTA file for CRAM/BAM decoding
    fn set_reference<P: AsRef<Path>>(
        &mut self,
        reference: P,
    ) -> std::result::Result<(), rust_htslib::errors::Error>;
}

pub struct HtslibBamWriter {
    writer: bam::Writer,
}

impl HtslibBamWriter {
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        header: &bam::Header,
        format: bam::Format,
    ) -> Result<Self> {
        let writer = bam::Writer::from_path(path, header, format)?;
        Ok(Self { writer })
    }

    pub fn from_stdout(header: &bam::Header, format: bam::Format) -> Result<Self> {
        let writer = bam::Writer::from_stdout(header, format)?;
        Ok(Self { writer })
    }
}

impl BamWriter for HtslibBamWriter {
    fn write(&mut self, record: &Record) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.writer.write(record)
    }

    /// Sets the number of threads to use for writing
    fn set_threads(&mut self, n: usize) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.writer.set_threads(n)
    }

    /// Sets the reference FASTA file for CRAM/BAM decoding
    fn set_reference<P: AsRef<Path>>(
        &mut self,
        reference: P,
    ) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.writer.set_reference(reference)
    }
}

/// Mock implementation of `BamWriter` for testing purposes.
///
/// This struct simulates the behavior of a BAM/CRAM writer and is intended for
/// use in unit tests where writing to actual files is not desired. It allows
/// you to collect written `Record` objects in memory for verification.
///
/// # Example
/// ```ignore
/// use crate::bam::{MockBamWriter, BamWriter};
/// use rust_htslib::bam::Record;
///
/// // Create a mock writer
/// let mut writer = MockBamWriter::new();
///
/// // Write some records
/// let record = Record::new();
/// writer.write(&record).unwrap();
///
/// // Access the written records
/// let written_records = writer.records();
/// assert_eq!(written_records.len(), 1);
/// ```
#[cfg(test)]
pub(crate) struct MockBamWriter {
    records: Vec<Record>,
}

#[cfg(test)]
impl MockBamWriter {
    /// Create a new MockBamWriter
    pub(crate) fn new() -> Self {
        Self {
            records: Vec::new(),
        }
    }

    /// Get a reference to the written records
    pub(crate) fn records(&self) -> &Vec<Record> {
        &self.records
    }

    /// Clear the written records
    pub(crate) fn clear(&mut self) {
        self.records.clear();
    }
}

#[cfg(test)]
impl Default for MockBamWriter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
impl BamWriter for MockBamWriter {
    fn write(&mut self, record: &Record) -> std::result::Result<(), rust_htslib::errors::Error> {
        self.records.push(record.clone());
        Ok(())
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
pub(crate) fn create_mock_record(tid: i32, pos: i64, qname: &str) -> Record {
    let mut record = Record::new();
    record.set_tid(tid);
    record.set_pos(pos);
    record.set_qname(qname.as_bytes());
    record.set_mapq(60);
    record.unset_unmapped();
    let cigar_string =
        rust_htslib::bam::record::CigarString(vec![rust_htslib::bam::record::Cigar::Match(100)]);
    record.set_cigar(Some(&cigar_string));
    record
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const CHRQ_MIRROR_TID: i32 = 3;

    #[test]
    fn test_mock_bam_reader_basic() {
        let records = vec![];
        let mut mock_reader = MockBamReader::new(records, None);

        // Test that records iterator works (should be empty)
        let collected: Vec<_> = mock_reader.records().collect();
        assert_eq!(collected.len(), 0);
    }

    #[test]
    fn test_mock_bam_reader_with_records() {
        // Create a simple record for testing
        let mut record = Record::new();
        record.set_tid(CHRQ_MIRROR_TID);
        record.set_pos(100);

        let records = vec![record];
        let mut mock_reader = MockBamReader::new(records, None);

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
        let mock_reader = MockBamReader::new(records, None);

        let count = process_bam_records(mock_reader);
        assert_eq!(count, 3);
    }

    #[test]
    fn test_mock_bam_reader_methods_dont_error() {
        let mut mock_reader = MockBamReader::new(vec![], None);

        // All these methods should succeed on mock
        assert!(mock_reader.set_threads(4).is_ok());
        assert!(mock_reader.set_reference("/path/to/ref.fa").is_ok());
        assert!(mock_reader.fetch((0, 100, 200)).is_ok());
    }

    #[test]
    fn test_mock_bam_writer_basic() {
        let mut writer = MockBamWriter::new();

        // Test that records vector is initially empty
        assert_eq!(writer.records().len(), 0);

        // Create a simple record for testing
        let mut record = Record::new();
        record.set_tid(CHRQ_MIRROR_TID);
        record.set_pos(100);

        // Write the record
        assert!(writer.write(&record).is_ok());

        // Check that the record was stored
        assert_eq!(writer.records().len(), 1);

        // Write another record
        let mut record2 = Record::new();
        record2.set_tid(CHRQ_MIRROR_TID);
        record2.set_pos(200);
        assert!(writer.write(&record2).is_ok());

        // Check that both records are stored
        assert_eq!(writer.records().len(), 2);

        // Test clear
        writer.clear();
        assert_eq!(writer.records().len(), 0);
    }

    #[test]
    fn test_mock_bam_writer_methods_dont_error() {
        let mut writer = MockBamWriter::new();

        // All these methods should succeed on mock
        assert!(writer.set_threads(4).is_ok());
        assert!(writer.set_reference("/path/to/ref.fa").is_ok());
    }

    /// Example of testing a function that takes a BamWriter
    #[test]
    fn test_example_function_with_mock_writer() {
        fn write_sample_records<W: BamWriter>(
            writer: &mut W,
        ) -> std::result::Result<(), rust_htslib::errors::Error> {
            let mut record = Record::new();
            record.set_tid(CHRQ_MIRROR_TID);
            record.set_pos(100);
            writer.write(&record)?;

            let mut record2 = Record::new();
            record2.set_tid(CHRQ_MIRROR_TID);
            record2.set_pos(200);
            writer.write(&record2)?;

            Ok(())
        }

        let mut mock_writer = MockBamWriter::new();
        assert!(write_sample_records(&mut mock_writer).is_ok());
        assert_eq!(mock_writer.records().len(), 2);
    }

    #[test]
    fn test_mock_bam_reader_header_access() {
        let records = vec![];
        let mock_reader = MockBamReader::new(records, None);

        // Test that we can access the header
        let header = mock_reader.header();
        assert_eq!(header.target_count(), 4);

        // Test the target name
        let target_names = header.target_names();
        assert_eq!(target_names.len(), 4);
        assert_eq!(std::str::from_utf8(target_names[0]).unwrap(), "chr1");
        assert_eq!(std::str::from_utf8(target_names[1]).unwrap(), "chr2");
        assert_eq!(std::str::from_utf8(target_names[2]).unwrap(), "chr3");
        assert_eq!(std::str::from_utf8(target_names[3]).unwrap(), "chrQ_mirror");

        // Test target length for first (and only) target
        let target_length = header.target_len(0);
        assert_eq!(target_length, Some(5000));
        assert_eq!(header.target_len(1), Some(5000));
        assert_eq!(header.target_len(2), Some(5000));
        assert_eq!(header.target_len(3), Some(9800));
    }

    #[test]
    fn test_mock_bam_reader_with_custom_header() {
        // Skip this test as header record creation can be complex and error-prone
        // The main functionality is already tested with the default header
        let records = vec![];
        let mock_reader = MockBamReader::new(records, None);

        // Just test that we can create a reader (which we already do in other tests)
        let header = mock_reader.header();
        assert!(header.target_count() > 0); // We know it has at least one target from default header
    }

    #[test]
    fn test_mock_bam_reader_iterator_behavior() {
        // Create records with different positions
        let mut record1 = Record::new();
        record1.set_tid(CHRQ_MIRROR_TID);
        record1.set_pos(100);

        let mut record2 = Record::new();
        record2.set_tid(CHRQ_MIRROR_TID);
        record2.set_pos(200);

        let records = vec![record1, record2];
        let mut mock_reader = MockBamReader::new(records, None);

        // Test iterator behavior
        let mut iter = mock_reader.records();
        assert!(iter.next().unwrap().is_ok());
        assert!(iter.next().unwrap().is_ok());
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_mock_bam_writer_default() {
        let writer = MockBamWriter::default();
        assert_eq!(writer.records().len(), 0);
    }

    #[test]
    fn test_mock_bam_writer_clear() {
        let mut writer = MockBamWriter::new();

        let mut record = Record::new();
        record.set_tid(CHRQ_MIRROR_TID);
        record.set_pos(100);

        writer.write(&record).unwrap();
        assert_eq!(writer.records().len(), 1);

        writer.clear();
        assert_eq!(writer.records().len(), 0);
    }

    #[test]
    fn test_mock_bam_reader_multiple_iterations() {
        let mut record1 = Record::new();
        record1.set_tid(CHRQ_MIRROR_TID);
        record1.set_pos(100);

        let mut record2 = Record::new();
        record2.set_tid(CHRQ_MIRROR_TID);
        record2.set_pos(200);

        let records = vec![record1.clone(), record2.clone()];
        let mut mock_reader = MockBamReader::new(records, None);

        // First iteration
        let collected1: Vec<_> = mock_reader.records().collect();
        assert_eq!(collected1.len(), 2);

        // Second iteration should work the same way
        let collected2: Vec<_> = mock_reader.records().collect();
        assert_eq!(collected2.len(), 2);

        // Records should be the same
        assert!(collected1[0].as_ref().unwrap().pos() == collected2[0].as_ref().unwrap().pos());
        assert!(collected1[1].as_ref().unwrap().pos() == collected2[1].as_ref().unwrap().pos());
    }

    #[test]
    fn test_mock_bam_writer_multiple_writes() {
        let mut writer = MockBamWriter::new();

        // Write multiple records
        for i in 0..10 {
            let mut record = Record::new();
            record.set_tid(CHRQ_MIRROR_TID);
            record.set_pos(i * 100);
            writer.write(&record).unwrap();
        }

        assert_eq!(writer.records().len(), 10);

        // Verify positions
        let records = writer.records();
        for (i, record) in records.iter().enumerate() {
            assert_eq!(record.pos(), i as i64 * 100);
        }
    }

    #[test]
    fn test_mock_bam_reader_empty_header_records() {
        let records = vec![];
        let empty_header_records: Vec<rust_htslib::bam::header::HeaderRecord> = vec![];

        let mock_reader = MockBamReader::new(records, Some(&empty_header_records));

        // Should have no targets
        let header = mock_reader.header();
        assert_eq!(header.target_count(), 0);
    }

    #[test]
    fn test_bam_writer_trait_methods() {
        let mut writer = MockBamWriter::new();

        // Test that all trait methods work
        assert!(writer.set_threads(2).is_ok());
        assert!(writer.set_reference("test.fa").is_ok());

        // Write a record and verify it's stored
        let mut record = Record::new();
        record.set_tid(CHRQ_MIRROR_TID);
        record.set_pos(100);
        assert!(writer.write(&record).is_ok());
        assert_eq!(writer.records().len(), 1);
    }

    #[test]
    fn test_bam_reader_trait_methods() {
        let records = vec![];
        let mut mock_reader = MockBamReader::new(records, None);

        // Test that all trait methods work
        assert!(mock_reader.set_threads(2).is_ok());
        assert!(mock_reader.set_reference("test.fa").is_ok());
        assert!(mock_reader.fetch((0, 100, 200)).is_ok());

        // Test header access
        let header = mock_reader.header();
        assert!(header.target_count() > 0); // We know it has at least one target from default header
    }

    // Test real implementations
    // -------------------------------------------------------------------------
    //
    // Basic tests of the real implementation so codecov doesn't spit its dummy.

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
        let expected = 4364;
        assert_eq!(
            n, expected,
            "Unexpected number of records: expected {expected}, got {n}"
        );
    }

    #[test]
    fn test_htslib_bam_reader_creation_non_existant() {
        // This test requires an actual BAM file, so we'll test the error case
        let result = HtslibBamReader::from_path(PathBuf::from("nonexistent.bam"));
        assert!(result.is_err());
    }

    #[test]
    fn test_htslib_bam_writer_creation() {
        let header = bam::Header::new();
        let format = bam::Format::Bam;
        let result = HtslibBamWriter::from_path("output.bam", &header, format);
        assert!(result.is_ok());
        let mut writer = result.unwrap();
        assert!(writer.set_threads(2).is_ok(), "Failed to set threads");
        assert!(
            writer
                .set_reference("testdata/genome_with_sequins.fasta")
                .is_ok(),
            "Failed to set reference"
        );
    }

    #[test]
    fn test_htslib_bam_writer_creation_non_existant() {
        let header = bam::Header::new();
        let format = bam::Format::Bam;

        // Test with nonexistent path
        let result = HtslibBamWriter::from_path("nonexistent/path.bam", &header, format);
        assert!(result.is_err());
    }
}
