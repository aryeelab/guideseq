#!/bin/bash
# This script downloads a full GUIDE-Seq dataset and performs runs the analysis pipeline.
# It should be run from the test directory.

cd large_test

# Create an output directory with a commit id hash suffix
OUTDIR=output.`git log --pretty=format:'%h' -n 1`
mkdir -p $OUTDIR
ln -sf $OUTDIR output

# Install bwa
git clone https://github.com/lh3/bwa.git
cd bwa
git checkout tags/0.7.9a
make
cd ..
PATH=`pwd`/bwa:$PATH

# Install bedtools
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
git checkout tags/v2.25.0
make
cd ..
PATH=`pwd`/bedtools2/bin:$PATH

# Download test data FASTQs and manifest
wget https://storage.googleapis.com/aryeelab/guideseq/guideseq_test_fastq.zip
unzip guideseq_test_fastq.zip

# Download the reference genome
wget http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta

# Run analysis pipeline
python ../../guideseq/guideseq.py all -m test_manifest.yaml

# Check that output tables match the reference output
cd output/filtered
md5sum -c ../../reference_output/md5.txt
