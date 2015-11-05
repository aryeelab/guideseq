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
TEST_DATA_FILES = {'forward': 'data/undemultiplexed/undemux.r1.fastq.gz',
                  'reverse': 'data/undemultiplexed/undemux.r2.fastq.gz',
                  'index1': 'data/undemultiplexed/undemux.i1.fastq.gz',
                  'index2': 'data/undemultiplexed/undemux.i2.fastq.gz'}
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

TEST_OUTPUT_PATH = 'test_output'
TEST_MIN_READS = 1000
TEST_DEMULTIPLEX_MANIFEST_PATH = os.path.join(TEST_OUTPUT_PATH, 'demultiplex_manifest.yaml')
TEST_MANIFEST_PATH = os.path.join(TEST_OUTPUT_PATH, 'test_manifest.yaml')

CORRECT_DEMULTIPLEXED_OUTPUT = 'data/demultiplexed'
CORRECT_UMITAGGED_OUTPUT = 'data/demultiplexed'
CORRECT_CONSOLDIATED_OUTPUT = 'data/demultiplexed'
CORRECT_ALIGNED_OUTPUT = 'data/demultiplexed'

class DemultiplexTestCase(unittest.TestCase):

    def setUp(self):
        # Create the test output folder
        os.makedirs(TEST_OUTPUT_PATH)

        # Create the test demultiplexing YAML
        test_manifest_data = {}
        test_manifest_data['undemultiplexed'] = TEST_DATA_FILES
        test_manifest_data['demultiplex_min_reads'] = TEST_MIN_READS
        test_manifest_data['samples'] = TEST_SAMPLES
        test_manifest_data['output_folder'] = TEST_OUTPUT_PATH

        with open(TEST_DEMULTIPLEX_MANIFEST_PATH, 'w') as f:
            f.write(yaml.dump(test_manifest_data, default_flow_style=False))


    def testIfDemultiplexed(self):
        g = guideseq.GuideSeq()
        g.parseManifestDemultiplex(TEST_DEMULTIPLEX_MANIFEST_PATH)
        g.demultiplex()

        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'demultiplexed'),
                                                  CORRECT_DEMULTIPLEXED_OUTPUT))

    def tearDown(self):
        # Delete temp output
        shutil.rmtree(TEST_OUTPUT_PATH)
        pass


class ConsolidateTestCase(unittest.TestCase):

    def setUp(self):
        # do the consolidation
        pass

    def testIfConsolidated(self):
        # Test if output equals consolidated output
        pass

    def tearDown(self):
        # Delete temp output
        pass


class AlignmentTestCase(unittest.TestCase):

    def setUp(self):
        # do the alignment
        pass

    def testIfAligned(self):
        # Test if output equals expected alignment output
        pass

    def tearDown(self):
        # Delete temp output
        pass


class OfftargetIdentificationTestCase(unittest.TestCase):

    def setUp(self):
        # do the offtarget identification
        pass

    def testIfConsolidated(self):
        # Test if output equals expected offtarget identification output
        pass

    def tearDown(self):
        # Delete temp output
        pass


class BackgroundSubtractTestCase(unittest.TestCase):

    def setUp(self):
        # do the bedtools subtraction
        pass

    def testIfConsolidated(self):
        # Test if output equals consolidated output
        pass

    def tearDown(self):
        # Delete temp output
        pass


if __name__ == '__main__':
    unittest.main()