import subprocess
import os


def filterBackgroundSites(bedtools_path, sample_name, sample_path, control_path, outfile):

    print 'Running background filtering for {0} sample'.format(sample_name)
    bedtools_filter_command = 'bedtools intersect -a {0} -b {1}'.format(sample_path, control_path)

    with open(outfile, 'w') as outfile:
        subprocess.call(bedtools_filter_command.split(), stdout=outfile)

    print 'Background filtering for {0} sample completed.'.format(sample_name)