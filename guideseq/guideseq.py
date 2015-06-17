# -*- coding: utf-8 -*-
"""

guideseq.py serves as the wrapper

"""

import yaml


class GuideSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path):

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.load(f)

        try:
            self.BWA_path  = manifest_data['BWA']
            self.HG19_path = manifest_data['HG19']
            self.samples   = manifest_data['samples']
        except:
            print 'Incomplete or incorrect manifest file. Please ensure your manifest contains all required fields.'

        if 'control' not in self.samples.keys():
            raise AssertionError('Your manifest must have a control sample specified.')

        if len(self.samples.keys()) < 2:
            raise AssertionError('Your manifest must have at least one control and one treatment sample.')


    def alignReads(self):
        import alignReads
        

    def demultiplex(self):
        pass

    def consolidate(self):
        pass

    def filter(self):
        pass

    def identifyOfftargetSites(self):
        pass

def main():
    g = GuideSeq()
    g.parseManifest('../tests/manifest.yaml')

    print g.samples

if __name__ == '__main__':
    main()