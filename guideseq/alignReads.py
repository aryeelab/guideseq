"""
alignReads
"""

import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False


def alignReads(BWA_path, HG19_path, read1, read2, outfile):

    sample_name = os.path.basename(outfile).split('.')[0]
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    sample_alignment_paths = {}

    # Check if genome is already indexed by bwa
    index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']

    genome_indexed = True
    for extension in index_files_extensions:
        if not os.path.isfile(HG19_path + extension):
            genome_indexed = False
            break

    # If the genome is not already indexed, index it
    if not genome_indexed:
        logger.info('Genome index files not detected. Running BWA to generate indices.')
        bwa_index_command = '{0} index {1}'.format(BWA_path, HG19_path)
        logger.info('Running bwa command: %s', bwa_index_command)
        subprocess.call(bwa_index_command.split())
        logger.info('BWA genome index generated')
    else:
        logger.info('BWA genome index found.')

    # Run paired end alignment against the genome
    logger.info('Running paired end mapping for {0}'.format(sample_name))
    bwa_alignment_command = '{0} mem {1} {2} {3}'.format(BWA_path,
                                                         HG19_path,
                                                         read1,
                                                         read2)

    logger.info(bwa_alignment_command)

    # Open the outfile and redirect the output of the alignment to it.
    with open(outfile, 'w') as f:
        subprocess.call(bwa_alignment_command.split(), stdout=f)

    logger.info('Paired end mapping for {0} completed.'.format(sample_name))
