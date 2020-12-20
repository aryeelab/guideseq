# This script generates a subsetted genome index for use in unit tests
# The index is hosted at: http://aryee.mgh.harvard.edu/guideseq/data/Homo_sapiens.GRCh38.dna.subset.masked.fa.index.zip
# Requirements: samtools, bedtools

mkdir -p genome_prep
cd genome_prep

# Download chromosomes 1 2 3 6 15
wget ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
wget ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
wget ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
wget ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz

cat *.fa.gz > Homo_sapiens.GRCh38.dna.subset.fa.gz
gunzip Homo_sapiens.GRCh38.dna.subset.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.subset.fa

# Pad test regions with 1kb on either side
bedtools slop -i ../test_regions.bed -g Homo_sapiens.GRCh38.dna.subset.fa.fai -b 1000 > test_regions_padded.bed

# Generate complement bed file (i.e. non-test regions)
bedtools complement -i test_regions_padded.bed -g Homo_sapiens.GRCh38.dna.subset.fa.fai > test_regions_complement.bed

# Mask non-test regions with Ns
bedtools maskfasta -fi Homo_sapiens.GRCh38.dna.subset.fa -fo Homo_sapiens.GRCh38.dna.subset.masked.fa -bed test_regions_complement.bed

# Move genome fasta to test dir
mv Homo_sapiens.GRCh38.dna.subset.masked.fa ..
cd ..

# Get bwa
git clone https://github.com/lh3/bwa.git
cd bwa
git checkout tags/0.7.12
make
cd ..

# Index the genome
bwa/bwa index Homo_sapiens.GRCh38.dna.subset.masked.fa