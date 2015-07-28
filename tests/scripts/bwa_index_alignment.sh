# Run bwa to generate index

bwa index ~/GRCh37/Homo_sapiens_assembly19.fasta

# Run paired end mapping to generate SAM files

bwa mem ~/GRCh37/Homo_sapiens_assembly19.fasta emx1.r1.fastq.gz emx1.r2.fastq.gz > ../output/emx1.sam
bwa mem ~/GRCh37/Homo_sapiens_assembly19.fasta control.r1.fastq.gz control.r2.fastq.gz > ../output/control.sam

