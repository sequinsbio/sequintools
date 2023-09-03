# Calibrating Sequin BAM Files

The original `sequintools` only has a very basic calibration tool that
downsamples a sequin region to a user specified coverage. The same target
coverage is used across the sequin region and is the same for all sequin regions
being processed. This tool also only outputs the sequin regions that are
calibrated, that is, any sample data is lost in the calibration process.

The problem with this approach to calibration is that it can result in coverage
that does not accurately reflect a samples actual coverage. This can be seen in the sequin
PMS2_ClinVar_1454389; this sequin contains the pathogenic variant
chr7:g.5973471_5973473delinsAACT and coverage the region chr7:5971971-5974970.
PMS2 has numerous pseudogenes and the detection of 3' variants is particularly
complicated by the pseudogene PMS2CL. NGS coverage of the 3' region of the
sequin region in GRCh38 is non-existant as everything mapping there has a mapQ
of 0. However, sequin data (on chrQ_mirror) does not suffer from this problem as
there is no PMS2CL corresponding region in the chrQ_mirror sequence all reads
generated from the sequin are aligned to their correct positions. This results
in uniform coverage across the sequin which is not an accurate reflection of the
real sample data.

The calibrate algorithm implemented in this program attempts to handle
situations like this by using the sample coverage to more accurately calibrate
the sequin coverage.


