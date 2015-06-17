import subprocess

def alignReads(BWA_path, HG19_path, samples):

    print 'Running BWA to generate index'
    subprocess.call('{0} index {1}'.format(BWA_path, HG19_path))
    print 'BWA genome index generated'

    for (sample_name, sample_paths) in samples.items():
        print 'Running paired end mapping for {0} sample'.format(sample_name)
        subprocess.call('{0} mem {1} {2} {3} > {4}'.format(BWA_path, HG19_path, 
                                                           sample_paths['forward'],
                                                           sample_paths['reverse'], 
                                                           os.path.join(out_path, sample_name + '.sam'
                                                           ))

        print 'Paired end mapping for {0} sample completed.'.format(sample_name)


"""
# Run bwa to generate index

bwa index ~/GRCh37/Homo_sapiens_assembly19.fasta

# Run paired end mapping to generate SAM files

bwa mem ~/GRCh37/Homo_sapiens_assembly19.fasta emx1.r1.fastq.gz emx1.r2.fastq.gz > ../output/emx1.sam
bwa mem ~/GRCh37/Homo_sapiens_assembly19.fasta control.r1.fastq.gz control.r2.fastq.gz > ../output/control.sam
"""
