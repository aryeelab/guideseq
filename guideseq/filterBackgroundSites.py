import subprocess
import os


def filterBackgroundSites(bedtools_path, samples, output_folder_path):

    sample_filtered_paths = {}

    for (sample_name, )

    # Run paired end alignment against the genome for each sample
    for (sample_name, sample_paths) in samples.items():
        print 'Running paired end mapping for {0} sample'.format(sample_name)
        bwa_alignment_command = '{0} mem {1} {2} {3}'.format(BWA_path, HG19_path,
                                                             sample_paths['forward'],
                                                             sample_paths['reverse'])

        # Open the outfile and redirect the output of the alignment to it.
        outfile_path = os.path.join(output_path, sample_name + '.sam')

        with open(outfile_path, 'w') as outfile:
            subprocess.call(bwa_alignment_command.split(), stdout=outfile)
            sample_alignment_paths[sample_name] = outfile_path

        print 'Paired end mapping for {0} sample completed.'.format(sample_name)

    return sample_alignment_paths
