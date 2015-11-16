#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_guideseq
----------------------------------

Tests for `guideseq` module.
"""

import yaml
import unittest
import os
import shutil
import utils
from guideseq import guideseq

TEST_SAMPLE_BARCODES = {'AGGCATGAGATCGC': 'mysample', 'GACTCCTGCGATAT': 'sample2'}
TEST_UNDEMULTIPLEXED_FILES = {'forward': 'data/undemultiplexed/undemux.r1.fastq.gz',
                  'reverse': 'data/undemultiplexed/undemux.r2.fastq.gz',
                  'index1': 'data/undemultiplexed/undemux.i1.fastq.gz',
                  'index2': 'data/undemultiplexed/undemux.i2.fastq.gz'}
TEST_DEMULTIPLEXED_FILES = {'read1': 'data/demultiplexed/EMX1.r1.fastq',
                            'read2': 'data/demultiplexed/EMX1.r2.fastq',
                            'index1': 'data/demultiplexed/EMX1.i1.fastq',
                            'index2': 'data/demultiplexed/EMX1.i2.fastq'}
TEST_SAMPLES = {
                'control':{
                 'barcode1':'CTCTCTAC',
                 'description':'Control',
                 'barcode2':'CTCTCTAT',
                 'target':None
                },
                'EMX1':{
                 'barcode1':'TAGGCATG',
                 'description':'EMX_site1',
                 'barcode2':'TAGATCGC',
                 'target':'GAGTCCGAGCAGAAGAAGAANGG'
                }
               }

TEST_SAMPLE_NAME = 'EMX1'
TEST_OUTPUT_PATH = 'test_output'
TEST_MIN_READS = 1000
TEST_DEMULTIPLEX_MANIFEST_PATH = os.path.join(TEST_OUTPUT_PATH, 'demultiplex_manifest.yaml')
TEST_MANIFEST_PATH = os.path.join(TEST_OUTPUT_PATH, 'test_manifest.yaml')

TEST_BWA_PATH = 'bwa'
TEST_BEDTOOLS_PATH = 'bedtools'

TEST_REFERENCE_GENOME = 'test_genome.fa'

CORRECT_DEMULTIPLEXED_OUTPUT = 'data/demultiplexed'
CORRECT_UMITAGGED_OUTPUT = 'data/umitagged'
CORRECT_CONSOLDIATED_OUTPUT = 'data/consolidated'
CORRECT_ALIGNED_OUTPUT = 'data/aligned'
CORRECT_IDENTIFIED_OUTPUT = 'data/identified'
CORRECT_FILTERED_OUTPUT = 'data/filtered'

CORRECT_ALL_OUTPUT = 'data'

class FullPipelineTestCase(unittest.TestCase):

    def setUp(self):
        # Create the test output folder
        os.makedirs(TEST_OUTPUT_PATH)

        # Create the test demultiplexing YAML
        test_manifest_data = {}
        test_manifest_data['undemultiplexed'] = TEST_UNDEMULTIPLEXED_FILES
        test_manifest_data['demultiplex_min_reads'] = TEST_MIN_READS
        test_manifest_data['samples'] = TEST_SAMPLES
        test_manifest_data['output_folder'] = TEST_OUTPUT_PATH
        test_manifest_data['bwa'] = TEST_BWA_PATH
        test_manifest_data['bedtools'] = TEST_BEDTOOLS_PATH
        test_manifest_data['reference_genome'] = TEST_REFERENCE_GENOME

        with open(TEST_MANIFEST_PATH, 'w') as f:
            f.write(yaml.dump(test_manifest_data, default_flow_style=False))


    def testFullPipeline(self):
        g = guideseq.GuideSeq()
        g.parseManifest(TEST_MANIFEST_PATH)

        # Demultiplex and test the demultiplex output
        g.demultiplex()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'demultiplexed'), CORRECT_DEMULTIPLEXED_OUTPUT))

        # UMITag and test the umitagging output
        g.umitag()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'umitagged'), CORRECT_UMITAGGED_OUTPUT))

        # Consolidate and test the consolidation output
        g.consolidate()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'consolidated'), CORRECT_CONSOLDIATED_OUTPUT))

        # Align and test the alignment output
        g.alignReads()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'aligned'), CORRECT_ALIGNED_OUTPUT))

        # Identify offtargets and test the output
        g.identifyOfftargetSites()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'identified'), CORRECT_IDENTIFIED_OUTPUT))

        # Filter background sites and test if correct
        g.filterBackgroundSites()
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'filtered'), CORRECT_FILTERED_OUTPUT))


    def tearDown(self):
        # Delete temp output
        shutil.rmtree(TEST_OUTPUT_PATH)
        pass

if __name__ == '__main__':
    unittest.main()