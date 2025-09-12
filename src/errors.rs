use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("BAM file not found")]
    BamFileNotFound,
    #[error("Invalid region specified")]
    InvalidRegion,
    #[error("An unknown error occurred: {0}")]
    Unknown(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("HTSlib error: {0}")]
    Hts(#[from] rust_htslib::errors::Error),

    #[error("UTF-8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),

    #[error("invalid BED record: {msg}")]
    BedInvalidRecord { msg: String },

    #[error("bedcov error: {msg}")]
    Bedcov { msg: String },

    #[error("calibration error: {msg}")]
    Calibration { msg: String },
}

pub type Result<T> = std::result::Result<T, Error>;
