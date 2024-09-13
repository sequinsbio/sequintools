# sequintools2

This is a reimplementation of the original `sequintools` in rust with
improvements.

## Installation

**(currently these don't actually work, but it would be cool if they did)**

If you're a Rust programmer, `sequintools` can be installed with `cargo`.

```sh
$ cargo install sequintools
```

Alternatively, one can use `cargo binstall` to install a `sequintools` binary
directly from GitHub:

```sh
$ cargo binstall sequintools
```

## Building

`sequintools` is written in Rust, so you'll need to grab a [Rust installation](https://www.rust-lang.org/) in order to compile it.

To build `sequintools`:

```sh
$ git clone https://github.com/sequinsbio/sequintools
$ cd sequintools
$ cargo build --release
$ sudo cp target/release/sequintools /usr/local/bin
```

## Usage 

Calibrate all sequin regions to a standard coverage of 40X

```sh
sequintools calibrate -b sequin_regions.chrQ_mirror.bed -f 40 -o calibrated.bam input.bam
```

Calibrate sequin coverage to mimic the coverage of the sample data and include
all sample data in output BAM

```sh
sequintools calibrate -b sequin_regions.chrQ_mirror.bed -H sequin_regions.hg38.bed \
    -o calibrated.bam input.bam
```

Same but exclude the sample data so that the output BAM only has calibrated
sequin data (this is much faster than the previous command).

```sh
sequintools calibrate -b sequin_regions.chrQ_mirror.bed -H sequin_regions.hg38.bed \
    -o calibrated.bam --exclude-uncalibrated-reads input.bam
```
