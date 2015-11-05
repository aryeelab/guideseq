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

with open('../.genome') as f:
    TEST_REFERENCE_GENOME = f.readline().strip()

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
        g.demultiplex()
        g.umitag()
        g.consolidate()
        g.alignReads()
        g.identifyOfftargetSites()
        g.filterBackgroundSites()

        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'demultiplexed'), CORRECT_DEMULTIPLEXED_OUTPUT))
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'umitagged'), CORRECT_UMITAGGED_OUTPUT))
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'consolidated'), CORRECT_CONSOLDIATED_OUTPUT))
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'aligned'), CORRECT_ALIGNED_OUTPUT))
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'identified'), CORRECT_IDENTIFIED_OUTPUT))
        self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'filtered'), CORRECT_FILTERED_OUTPUT))

    def tearDown(self):
        # Delete temp output
        shutil.rmtree(TEST_OUTPUT_PATH)
        pass

# class DemultiplexTestCase(unittest.TestCase):
#
#     def setUp(self):
#         # Create the test output folder
#         os.makedirs(TEST_OUTPUT_PATH)
#
#         # Create the test demultiplexing YAML
#         test_manifest_data = {}
#         test_manifest_data['undemultiplexed'] = TEST_UNDEMULTIPLEXED_FILES
#         test_manifest_data['demultiplex_min_reads'] = TEST_MIN_READS
#         test_manifest_data['samples'] = TEST_SAMPLES
#         test_manifest_data['output_folder'] = TEST_OUTPUT_PATH
#
#         with open(TEST_DEMULTIPLEX_MANIFEST_PATH, 'w') as f:
#             f.write(yaml.dump(test_manifest_data, default_flow_style=False))
#
#
#     def testIfDemultiplexed(self):
#         g = guideseq.GuideSeq()
#         g.parseManifestDemultiplex(TEST_DEMULTIPLEX_MANIFEST_PATH)
#         g.demultiplex()
#
#         self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'demultiplexed'),
#                                                   CORRECT_DEMULTIPLEXED_OUTPUT))
#
#     def tearDown(self):
#         # Delete temp output
#         shutil.rmtree(TEST_OUTPUT_PATH)


# class UMITagTestCase(unittest.TestCase):
#
#     def setUp(self):
#         # Create the test output folder
#         os.makedirs(TEST_OUTPUT_PATH)
#
#     def testIfUmiTagged(self):
#         # Test if output equals consolidated output
#         g = guideseq.GuideSeq()
#         g.output_folder = TEST_OUTPUT_PATH
#         g.samples = [TEST_SAMPLE_NAME]
#         g.demultiplexed = {TEST_SAMPLE_NAME: {}}
#         g.demultiplexed[TEST_SAMPLE_NAME] = TEST_DEMULTIPLEXED_FILES
#         g.umitag()
#
#         self.assertTrue(utils.checkFolderEquality(os.path.join(TEST_OUTPUT_PATH, 'umitagged'), CORRECT_UMITAGGED_OUTPUT))
#
#     def tearDown(self):
#         # Delete temp output
#         shutil.rmtree(TEST_OUTPUT_PATH)
#
# class ConsolidateTestCase(unittest.TestCase):
#
#     def setUp(self):
#         # Create the test output folder
#         os.makedirs(TEST_OUTPUT_PATH)
#
#     def testIfConsolidated(self):
#         # Test if output equals consolidated output
#
#         pass
#
#     def tearDown(self):
#         # Delete temp output
#         shutil.rmtree(TEST_OUTPUT_PATH)
#
#
# class AlignmentTestCase(unittest.TestCase):
#
#     def setUp(self):
#         # do the alignment
#         pass
#
#     def testIfAligned(self):
#         # Test if output equals expected alignment output
#         pass
#
#     def tearDown(self):
#         # Delete temp output
#         pass
#
#
# class OfftargetIdentificationTestCase(unittest.TestCase):
#
#     def setUp(self):
#         # do the offtarget identification
#         pass
#
#     def testIfConsolidated(self):
#         # Test if output equals expected offtarget identification output
#         pass
#
#     def tearDown(self):
#         # Delete temp output
#         pass
#
#
# class BackgroundSubtractTestCase(unittest.TestCase):
#
#     def setUp(self):
#         # do the bedtools subtraction
#         pass
#
#     def testIfConsolidated(self):
#         # Test if output equals consolidated output
#         pass
#
#     def tearDown(self):
#         # Delete temp output
#         pass
#

if __name__ == '__main__':
    unittest.main()