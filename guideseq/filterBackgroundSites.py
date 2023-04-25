import subprocess
import os

def filterBackgroundSites(bedtools_path, sample_path, control_path, outfile):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    sample_noHeader = os.path.join(os.path.dirname(sample_path), 'sample_noHeader.txt')
    control_noHeader = os.path.join(os.path.dirname(control_path), 'control_noHeader.txt')

    sample_noHeader_command = "sed '1d' {0} > {1}".format(sample_path, sample_noHeader)
    control_noHeader_command = "sed '1d' {0} > {1}".format(control_path, control_noHeader)
    clean_command = "rm {0} {1}".format(control_noHeader, sample_noHeader)
    bedtools_filter_command = '{0} intersect -a {1} -b {2}'.format(bedtools_path, sample_noHeader, control_noHeader)

    subprocess.check_call(sample_noHeader_command, shell=True, env=os.environ.copy())
    subprocess.check_call(control_noHeader_command, shell=True, env=os.environ.copy())

    with open(outfile, 'w') as output_file:
        subprocess.check_call(bedtools_filter_command, shell=True, env=os.environ.copy(), stdout=output_file)
    subprocess.check_call(clean_command, shell=True, env=os.environ.copy())
