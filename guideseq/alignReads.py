import subprocess
import os


def alignReads(BWA_path, HG19_path, samples, output_path):

    # Check if genome is already indexed by bwa
    index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']

    genome_indexed = True
    for extension in index_files_extensions:
        if not os.path.isfile(HG19_path + extension):
            genome_indexed = False
            break

    # If the genome is not already indexed, index it
    if not genome_indexed:
        print 'Genome index files not detected. Running BWA to generate indices.'
        bwa_index_command = '{0} index {1}'.format(BWA_path, HG19_path)
        subprocess.call(bwa_index_command.split())
        print 'BWA genome index generated'

    # Run paired end alignment against the genome for each sample
    for (sample_name, sample_paths) in samples.items():
        print 'Running paired end mapping for {0} sample'.format(sample_name)
        bwa_alignment_command = '{0} mem {1} {2} {3}'.format(BWA_path, HG19_path,
                                                            sample_paths['forward'],
                                                            sample_paths['reverse'])

        # Open the outfile and redirect the output of the alignment to it.
        outfile_path = os.path.join(output_path, sample_name + '.sam')
        with open(outfile_path) as outfile:
            subprocess.call(bwa_alignment_command.split(), stdout=outfile)

        print 'Paired end mapping for {0} sample completed.'.format(sample_name)
