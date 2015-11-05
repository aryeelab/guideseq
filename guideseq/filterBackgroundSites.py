import subprocess
import os


def filterBackgroundSites(bedtools_path, sample_path, control_path, outfile):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    bedtools_filter_command = 'bedtools intersect -a {0} -b {1}'.format(sample_path, control_path)

    with open(outfile, 'w') as outfile:
        subprocess.call(bedtools_filter_command.split(), stdout=outfile)
