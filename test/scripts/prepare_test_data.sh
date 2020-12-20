## This script generates two test FASTQ datasets:
##  1. An undemultiplexed dataset representing several barcoded samples with molecular indexing.
##     This dataset represents raw data from a MiSeq run following the GUIDE-Seq protocol
##     described in Tsai et al., 2014 (PMID XXXX)
##  2. A two sample dataset containing reads from a control and an EMX guide experiment
##     that overlap with a small set of test regions: 
##     the on-target location, 3 off-target locations and two DSB hotspots. 
##     The reads representing the same template molecule (i.e. those with the same 
##     molecular barcode have been consolidated).
##
## Raw input dataset:
## /data/joung/sequencing_fastq/131007_M01326_0075_000000000-A6B33/fastq_with_indexes
## -rw-rw----. 1 ma695 aryee 2.7G Oct 31 12:44 guideseq_test_fastq.zip
## -rw-rw----. 1 st680 joung 120M Oct 14 14:52 Undetermined_S0_L001_I1_001.fastq.gz
## -rw-rw----. 1 st680 joung 221M Oct 14 14:53 Undetermined_S0_L001_I2_001.fastq.gz
## -rw-rw----. 1 st680 joung 1.1G Oct 14 14:58 Undetermined_S0_L001_R1_001.fastq.gz
## -rw-rw----. 1 st680 joung 1.3G Oct 14 14:59 Undetermined_S0_L001_R2_001.fastq.gz
##
## EMX1 has barcode P706 (TAGGCATG), A01 (TAGATCGC).
## Ignoring the first base and concatenating gives AGGCATGAGATCGC.
## EMX1 target sequence: GAGTCCGAGCAGAAGAAGAANGG
##
## Oligo control has barcode P707 (CTCTCTAC), A02 (CTCTCTAT)
## Ignoring the first base and concatenating gives TCTCTACTCTCTAT.

ON_TARGET="2:73160981-73161004" 
OFF_TARGET="15:44109746-44109769 6:9118792-9118815 2:218378101-218378124"
DSB_HOTSPOTS="1:236260170-236260754 3:197900267-197900348"

UMI_PKG_DIR="../../guideseq/umi"

# Align reads
INPUT_DIR="/data/joung/sequencing_fastq/131007_M01326_0075_000000000-A6B33/fastq_with_indexes"
BWA_INDEX="/data/aryee/pub/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"
module load aryee/bwa-0.7.9a 
time bwa mem $BWA_INDEX $INPUT_DIR/Undetermined_S0_L001_R1_001.fastq.gz $INPUT_DIR/Undetermined_S0_L001_R2_001.fastq.gz > undemux.sam

# Generate BAM:
module load samtools/0.1.19
samtools view -bS undemux.sam > undemux.bam

# Sort BAM
samtools sort undemux.bam undemux.sorted

# Index BAMs
samtools index undemux.sorted.bam

# Get the names of reads that overlap with the selected test regions:
samtools view undemux.sorted.bam $ON_TARGET $OFF_TARGET $DSB_HOTSPOTS | cut -f1 | sort | uniq > read_names.txt

# Subset FASTQs to extract _all_ read pairs where at least one of the reads falls in a specified test region
zcat $INPUT_DIR/Undetermined_S0_L001_R1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > undemux_all.r1.fastq.gz
zcat $INPUT_DIR/Undetermined_S0_L001_R2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > undemux_all.r2.fastq.gz
zcat $INPUT_DIR/Undetermined_S0_L001_I1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > undemux_all.i1.fastq.gz
zcat $INPUT_DIR/Undetermined_S0_L001_I2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names.txt | gzip -c > undemux_all.i2.fastq.gz

# Demultiplex full target region FASTQs
python $UMI_PKG_DIR/demultiplex.py --min_reads 1000 --read1 undemux_all.r1.fastq.gz --read2 undemux_all.r2.fastq.gz --index1 undemux_all.i1.fastq.gz --index2 undemux_all.i2.fastq.gz --sample_barcodes samplekey.txt

# Choose a subset of EMX1 and control read names:
cat emx1.r1.fastq  | grep ^@M01326 | cut -f1 -d ' ' | sort | uniq | shuf --random-source emx1.r1.fastq -n 6000 > read_names_sample.txt
cat control.r1.fastq  | grep ^@M01326 | cut -f1 -d ' ' | sort | uniq | shuf --random-source control.r1.fastq -n 2000 >> read_names_sample.txt

# Subset FASTQs to extract _a sample of_ read pairs where at least one of the reads falls in a specified test region
zcat $INPUT_DIR/Undetermined_S0_L001_R1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names_sample.txt | gzip -c > undemux.r1.fastq.gz
zcat $INPUT_DIR/Undetermined_S0_L001_R2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names_sample.txt | gzip -c > undemux.r2.fastq.gz
zcat $INPUT_DIR/Undetermined_S0_L001_I1_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names_sample.txt | gzip -c > undemux.i1.fastq.gz
zcat $INPUT_DIR/Undetermined_S0_L001_I2_001.fastq.gz | grep -F -A3 --no-group-separator -f read_names_sample.txt | gzip -c > undemux.i2.fastq.gz

# Demultiplex sub-sampled target region FASTQs
python $UMI_PKG_DIR/demultiplex.py --min_reads 1000 --read1 undemux.r1.fastq.gz --read2 undemux.r2.fastq.gz --index1 undemux.i1.fastq.gz --index2 undemux.i2.fastq.gz --sample_barcodes samplekey.txt

# Consolidate reads with the same molecular index
for SAMPLE in emx1 control
do
    echo "Consolidating reads for $SAMPLE"
    python $UMI_PKG_DIR/umitag.py --read1_in $SAMPLE.r1.fastq --read2_in $SAMPLE.r2.fastq --read1_out $SAMPLE.r1.umitagged.fastq --read2_out $SAMPLE.r2.umitagged.fastq --index1 $SAMPLE.i1.fastq --index2 $SAMPLE.i2.fastq
    python $UMI_PKG_DIR/consolidate.py $SAMPLE.r1.umitagged.fastq $SAMPLE.r1.consolidated.fastq 15 0.9
    python $UMI_PKG_DIR/consolidate.py $SAMPLE.r2.umitagged.fastq $SAMPLE.r2.consolidated.fastq 15 0.9
done

# Copy test datasets to data dir
cp undemux.*.fastq.gz data
for SAMPLE in emx1 control
do
    gzip -c $SAMPLE.r1.consolidated.fastq > data/$SAMPLE.r1.fastq.gz
    gzip -c $SAMPLE.r2.consolidated.fastq > data/$SAMPLE.r2.fastq.gz
done
