# `sequintools`

`sequintools` is a suite of tools for manipulating and reporting on Next
Generation Sequencing data that has Sequins added to the sample.

## Installation

If you're a Rust programmer, `sequintools` can be installed with `cargo`.

```sh
cargo install sequintools
```

Alternatively, one can use `cargo binstall` to install a `sequintools` binary
directly from GitHub:

```sh
cargo binstall sequintools
```

**Note** that `sequintools` has a dependency on `cmake` and `clang`, and these
will need to be installed on your device before you can install `sequintools`.
On a Linux device you would usually install these using `apt`, e.g.

```sh
sudo apt update && sudo apt install -y cmake clang
```

If you're not using a Linux device then you'll need to research how to install
the dependencies on your device. Alternatively you can use the docker build
which will already have the dependencies required installed.

## Docker

We also provide `sequintools` as a [Docker container](https://github.com/orgs/sequinsbio/packages/container/package/sequintools).

This command line below is used to fetch its latest version.

```sh
docker pull ghcr.io/sequinsbio/sequintools:latest
```

As an example, you could run the calibration step as below, make sure the
example folder is under you current path.

```sh
docker run --rm -v $PWD:/data ghcr.io/sequinsbio/sequintools \
    sequintools calibrate \
    -b /data/example/resources/sequin_regions.chrQ_mirror.bed \
    -f 40 \
    -o /data/calibrated.bam \
    --write-index \
    /data/example/example.bam
```

## Building

`sequintools` is written in Rust, so you'll need to grab a [Rust
installation](https://www.rust-lang.org/) in order to compile it.

To build `sequintools`:

```sh
git clone https://github.com/sequinsbio/sequintools
cd sequintools
cargo build --release
sudo cp target/release/sequintools /usr/local/bin
```

## Usage

To get going, we provide some example data in the `example` directory. This is
simulated data that will allow you to quickly run commands from `sequintools`
and familiarise yourself with the input file requirements.

### `calibrate`

`calibrate` downsamples Sequin reads so that they more closely match the sample
data. The more closely the Sequins data resembles the sample data it controls
for, the better the control will be.

The simplest form of calibration is to adjust the mean coverage for all Sequin
regions to a standard coverage. In this example, the mean coverage for each
Sequin in the data set is set to 40X.

```sh
sequintools calibrate \
    -b example/resources/sequin_regions.chrQ_mirror.bed \
    -f 40 \
    -o calibrated.bam \
    --write-index \
    example/example.bam
```

Alternatively, you can use the sample data in the same BAM file to adjust the
Sequin coverage to more closely represent the coverage of the controlled
region. This method uses the mean depth of the region in the sample data that
corresponds to each Sequin rather than a fixed coverage for all Sequin regions.
To do this, you need to provide two BED files one containing the Sequin regions
on the decoy chromosome and the other containing the reference genome locations
they control for. The name column (4th column) must match between the two BED
files.

```sh
sequintools calibrate \
    -b example/resources/sequin_regions.chrQ_mirror.bed \
    -H example/resources/sequin_regions.hg38.bed \
    -o calibrated.bam \
    --write-index \
    example/example.bam
```

Make the same adjustments as the previous command, but exclude the sample data
so that the output BAM only has calibrated Sequin data (this is much faster than
the previous command).

```sh
sequintools calibrate \
    -b example/resources/sequin_regions.chrQ_mirror.bed \
    -H example/resources/sequin_regions.hg38.bed \
    -o calibrated.bam \
    --write-index \
    --exclude-uncalibrated-reads \
    example/example.bam
```

You can also create a summary file with the before and after coverage with the
`--summary-report` option:

```sh
sequintools calibrate \
    -b example/resources/sequin_regions.chrQ_mirror.bed \
    -H example/resources/sequin_regions.hg38.bed \
    -o calibrated.bam \
    --summary-report calibrate.summary.csv \
    --write-index \
    --exclude-uncalibrated-reads \
    example/example.bam
```

### `bedcov`

The `bedcov` command collects statistics from a BAM file for regions in a BED
file and writes them as a CSV. Statistics include the minimum and maximum
coverage per base, mean coverage, standard deviation of coverage and the
coefficient of variance.

```sh
sequintools bedcov \
    example/resources/sequin_regions.hg38.bed \
    example/example.bam
```
