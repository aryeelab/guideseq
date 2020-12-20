#!/bin/bash
# This script generates a subsetted genome index for use in unit tests
# Requirements: samtools, bedtools

mkdir -p genome_prep
cd genome_prep

# Download chromosomes 1 2 3 6 15
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.2.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.3.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.6.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.15.fa.gz

cat *.fa.gz > Homo_sapiens.GRCh37.75.dna.subset.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.subset.fa.gz
samtools faidx Homo_sapiens.GRCh37.75.dna.subset.fa

# Pad test regions with 1kb on either side
bedtools slop -i ../test_regions.bed -g Homo_sapiens.GRCh37.75.dna.subset.fa.fai -b 1000 > test_regions_padded.bed

# Extract test genome regions
bedtools getfasta -fi Homo_sapiens.GRCh37.75.dna.subset.fa -bed test_regions_padded.bed -fo test_genome.fa

# Move genome fasta to test dir
mv test_genome.fa ../..
cd ..

