#!/usr/bin/env bash
# This script is how the example simulated data was created.
set -euo pipefail

REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta

# Download and prepare the Human GRCh38 reference genome, if it doesn't already
# exist.
if ! [[ -f "$REF" ]]; then
	curl -fL https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz |
		gzip -dc >"$REF"
fi

# Create some simulated data based on the sequin set. Here we are simulating
# the sample data at ~30X and the sequin data at ~500X. This should give us
# plenty of data to downsample to our samples coverage.
art_illumina -i resources/sequin_sequences-natural.fa -ss HSXn -p -na -l 150 -f 30 -m 500 -s 20 --rndSeed 5678 -o nat_R
art_illumina -i resources/sequin_sequences-mirror.fa -ss HSXn -p -na -l 150 -f 100 -m 500 -s 20 --rndSeed 5678 -o sqn_R
cat nat_R1.fq sqn_R1.fq | gzip -c >sim_R1.fq.gz
cat nat_R2.fq sqn_R2.fq | gzip -c >sim_R2.fq.gz
rm nat_R1.fq nat_R2.fq sqn_R1.fq sqn_R2.fq

# Build a reference genome with GRCh38 plus the decoy sequence for the sequin set.
cat "$REF" resources/sequin_decoy.chrQ_mirror.fa >genome.fa
samtools faidx genome.fa
bwa index genome.fa

# Finally align our simulated data to the reference
bwa mem -M -t 8 -R '@RG\tID:A\tSM:A\tLB:A\tPU:A\tPL:ILLUMINA' genome.fa sim_R1.fq.gz sim_R2.fq.gz |
	samtools sort -o example.bam
samtools index example.bam

# You can view the results in IGV
# igv --genome genome.fa example.bam resources/sequin_regions.chrQ_mirror.bed resources/sequin_regions.hg38.bed
