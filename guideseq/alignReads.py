import subprocess

def alignReads(BWA_path, HG19_path, samples, output_path):

    print 'Running BWA to generate index'
    print '{0} index {1}'.format(BWA_path, HG19_path)
    subprocess.call('{0} index {1}'.format(BWA_path, HG19_path))
    print 'BWA genome index generated'

    for (sample_name, sample_paths) in samples.items():
        print sample_name
        print 'Running paired end mapping for {0} sample'.format(sample_name)
        subprocess.call('{0} mem {1} {2} {3} > {4}'.format(BWA_path, HG19_path,
                                                           sample_paths['forward'],
                                                           sample_paths['reverse'],
                                                           os.path.join(output_path, sample_name + '.sam')
                                                           ))

        print 'Paired end mapping for {0} sample completed.'.format(sample_name)