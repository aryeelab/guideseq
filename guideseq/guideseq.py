#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========
serves as the wrapper for all guideseq pipeline

"""

import os
import sys
import yaml
import argparse
import traceback

# Set up logger
import log
logger = log.createCustomLogger('root')

from alignReads import alignReads
from filterBackgroundSites import filterBackgroundSites
from umi import demultiplex, umitag, consolidate
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation

DEFAULT_DEMULTIPLEX_MIN_READS = 10000
DEFAULT_WINDOW_SIZE = 25
DEFAULT_MAX_SCORE = 7

CONSOLIDATE_MIN_QUAL = 15
CONSOLIDATE_MIN_FREQ = 0.9


class GuideSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path):
        logger.info('Loading manifest...')

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.safe_load(f)
        
        if not "cores" in manifest_data:
            manifest_data['cores'] = 4
        
        # Set default tag/primer sequences if not specified
        if not "primer1" in manifest_data:
            manifest_data['primer1'] = 'TTGAGTTGTCATATGTTAAT'
        if not "primer2" in manifest_data:
            manifest_data['primer2'] = 'ACATATGACAACTCAATTAA'

        try:
            # Validate manifest data
            validation.validateManifest(manifest_data)

            self.cores = manifest_data['cores']
            self.BWA_path = manifest_data['bwa']
            self.bedtools = manifest_data['bedtools']
            self.reference_genome = manifest_data['reference_genome']
            self.output_folder = manifest_data['output_folder']
            self.undemultiplexed = manifest_data['undemultiplexed']
            self.samples = manifest_data['samples']
            self.primer1 = manifest_data['primer1']
            self.primer2 = manifest_data['primer2']

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            sys.exit()

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in manifest_data:
            self.demultiplex_min_reads = manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS
        # Allow the user to specify window size for off-target search
        if 'search_radius' in manifest_data:
            self.search_radius = manifest_data['search_radius']
        else:
            self.search_radius = DEFAULT_WINDOW_SIZE
        # Allow the user to specify window size for off-target search
        if 'max_score' in manifest_data:
            self.max_score = manifest_data['max_score']
        else:
            self.max_score = DEFAULT_MAX_SCORE
        # Allow the user to specify PAM seq. Yichao 3/6/2020
        if 'PAM' in manifest_data:
            self.PAM = manifest_data['PAM']
        else:
            self.PAM = "NGG"

        # Make sure the user has specified a control barcode
        if 'control' not in self.samples.keys():
            raise AssertionError('Your manifest must have a control sample specified.')

        # Make sure the user has both a sample and a control
        if len(self.samples) < 2:
            raise AssertionError('Your manifest must have at least one control and one treatment sample.')

        logger.info('Successfully loaded manifest.')

    def parseManifestDemultiplex(self, manifest_path):
        logger.info('Loading manifest for demultiplexing...')

        with open(manifest_path, 'r') as f:
            manifest_data = yaml.load(f)

            try:
                self.output_folder = manifest_data['output_folder']
                self.undemultiplexed = manifest_data['undemultiplexed']
                self.samples = manifest_data['samples']

            except Exception as e:
                logger.error('Incomplete or incorrect manifest file. Please ensure your manifest contains all required fields.')
                quit()

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in manifest_data:
            self.demultiplex_min_reads = manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS

        logger.info('Successfully loaded manifest for single-step demultiplexing.')

    def demultiplex(self):

        logger.info('Demultiplexing undemultiplexed files...')

        # Take our two barcodes and concatenate them
        swapped_sample_barcodes = {}
        for sample in self.samples:
            barcode1 = self.samples[sample]['barcode1']
            barcode2 = self.samples[sample]['barcode2']
            barcode = barcode1[1:8] + barcode2[1:8]
            swapped_sample_barcodes[barcode] = sample

        try:
            demultiplex.demultiplex(self.undemultiplexed['forward'],
                                    self.undemultiplexed['reverse'],
                                    self.undemultiplexed['index1'],
                                    self.undemultiplexed['index2'],
                                    swapped_sample_barcodes,
                                    os.path.join(self.output_folder, 'demultiplexed'),
                                    min_reads=self.demultiplex_min_reads)

            self.demultiplexed = {}
            for sample in self.samples:
                self.demultiplexed[sample] = {}
                self.demultiplexed[sample]['read1'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.r1.fastq')
                self.demultiplexed[sample]['read2'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.r2.fastq')
                self.demultiplexed[sample]['index1'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.i1.fastq')
                self.demultiplexed[sample]['index2'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.i2.fastq')

            logger.info('Successfully demultiplexed reads.')
        except Exception as e:
            logger.error('Error demultiplexing reads.')
            logger.error(traceback.format_exc())
            quit()

    def umitag(self):
        logger.info('umitagging reads...')

        try:
            self.umitagged = {}
            for sample in self.samples:
                self.umitagged[sample] = {}
                self.umitagged[sample]['read1'] = os.path.join(self.output_folder, 'umitagged', sample + '.r1.umitagged.fastq')
                self.umitagged[sample]['read2'] = os.path.join(self.output_folder, 'umitagged', sample + '.r2.umitagged.fastq')

                umitag.umitag(self.demultiplexed[sample]['read1'],
                              self.demultiplexed[sample]['read2'],
                              self.demultiplexed[sample]['index1'],
                              self.demultiplexed[sample]['index2'],
                              self.umitagged[sample]['read1'],
                              self.umitagged[sample]['read2'],
                              os.path.join(self.output_folder, 'umitagged'))

            logger.info('Successfully umitagged reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def consolidate(self, min_freq=CONSOLIDATE_MIN_FREQ, min_qual=CONSOLIDATE_MIN_QUAL):
        logger.info('Consolidating reads...')

        try:
            self.consolidated = {}

            for sample in self.samples:
                self.consolidated[sample] = {}
                self.consolidated[sample]['read1'] = os.path.join(self.output_folder, 'consolidated', sample + '.r1.consolidated.fastq')
                self.consolidated[sample]['read2'] = os.path.join(self.output_folder, 'consolidated', sample + '.r2.consolidated.fastq')

                consolidate.consolidate(self.umitagged[sample]['read1'], self.consolidated[sample]['read1'], min_qual, min_freq)
                consolidate.consolidate(self.umitagged[sample]['read2'], self.consolidated[sample]['read2'], min_qual, min_freq)

            logger.info('Successfully consolidated reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def alignReads(self):
        logger.info('Aligning reads...')

        try:
            self.aligned = {}
            for sample in self.samples:
                sample_alignment_path = os.path.join(self.output_folder, 'aligned', sample + '.sam')
                alignReads(self.cores,
                           self.BWA_path,
                           self.reference_genome,
                           self.consolidated[sample]['read1'],
                           self.consolidated[sample]['read2'],
                           sample_alignment_path)
                self.aligned[sample] = sample_alignment_path
                logger.info('Finished aligning reads to genome.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()

    def identifyOfftargetSites(self):
        logger.info('Identifying offtarget sites...')

        try:
            self.identified = {}

            # Identify offtarget sites for each sample
            for sample in self.samples:

                # Prepare sample annotations
                sample_data = self.samples[sample]
                annotations = {}
                annotations['Description'] = sample_data['description']
                annotations['Targetsite'] = sample

                if sample == 'control':
                    annotations['Sequence'] = ''
                else:
                    annotations['Sequence'] = sample_data['target']

                samfile = os.path.join(self.output_folder, 'aligned', sample + '.sam')

                self.identified[sample] = os.path.join(self.output_folder, 'identified', sample + '_identifiedOfftargets.txt')

                identifyOfftargetSites.analyze(samfile, self.reference_genome, self.identified[sample], annotations,
                                               self.search_radius, self.max_score, self.primer1, self.primer2)

            logger.info('Finished identifying offtarget sites.')

        except Exception as e:
            logger.error('Error identifying offtarget sites.')
            logger.error(traceback.format_exc())
            quit()

    def filterBackgroundSites(self):
        logger.info('Filtering background sites')

        try:
            self.filtered = {}

            # Filter background in each sample
            for sample in self.samples:
                if sample != 'control':
                    self.filtered[sample] = os.path.join(self.output_folder, 'filtered', sample + '_backgroundFiltered.txt')
                    filterBackgroundSites(self.bedtools, self.identified[sample], self.identified['control'], self.filtered[sample])
                    logger.info('Finished background filtering for {0} sample'.format(sample))

            logger.info('Finished filtering background sites.')

        except Exception as e:
            logger.error('Error filtering background sites.')
            logger.error(traceback.format_exc())

    def visualize(self):
        logger.info('Visualizing off-target sites')

        # try:
            # for sample in self.samples:
                # if sample != 'control':
                    # infile = self.identified[sample]
                    # outfile = os.path.join(self.output_folder, 'visualization', sample + '_offtargets')
                    # visualizeOfftargets(infile, outfile, title=sample)

            # logger.info('Finished visualizing off-target sites')

        # except Exception as e:
            # logger.error('Error visualizing off-target sites.')
            # logger.error(traceback.format_exc())

        for sample in self.samples: ## 3/6/2020 Yichao solved: visualization stopped when one sample failed
            if sample != 'control':
                try:
                    infile = self.identified[sample]
                    outfile = os.path.join(self.output_folder, 'visualization', sample + '_offtargets')
                    try:
                        self.PAM
                        visualizeOfftargets(infile, outfile, title=sample,PAM=self.PAM)
                    except:
                        visualizeOfftargets(infile, outfile, title=sample,PAM="NGG")
                except Exception as e:
                    logger.error('Error visualizing off-target sites: %s'%(sample))
                    logger.error(traceback.format_exc())
        logger.info('Finished visualizing off-target sites')

def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    all_parser.add_argument('--identifyAndFilter', action='store_true', default=False)
    all_parser.add_argument('--skip_demultiplex', action='store_true', default=False)

    demultiplex_parser = subparsers.add_parser('demultiplex', help='Demultiplex undemultiplexed FASTQ files')
    demultiplex_parser.add_argument('--manifest', '-m', help='Specify the manifest path', required=True)

    umitag_parser = subparsers.add_parser('umitag', help='UMI tag demultiplexed FASTQ files for consolidation')
    umitag_parser.add_argument('--read1', required=True)
    umitag_parser.add_argument('--read2', required=True)
    umitag_parser.add_argument('--index1', required=True)
    umitag_parser.add_argument('--index2', required=True)
    umitag_parser.add_argument('--outfolder', required=True)

    consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate UMI tagged FASTQs')
    consolidate_parser.add_argument('--read1', required=True)
    consolidate_parser.add_argument('--read2', required=True)
    consolidate_parser.add_argument('--outfolder', required=True)
    consolidate_parser.add_argument('--min_quality', required=False, type=float)
    consolidate_parser.add_argument('--min_frequency', required=False, type=float)

    align_parser = subparsers.add_parser('align', help='Paired end read mapping to genome')
    align_parser.add_argument('--bwa', required=True)
    align_parser.add_argument('--genome', required=True)
    align_parser.add_argument('--read1', required=True)
    align_parser.add_argument('--read2', required=True)
    align_parser.add_argument('--outfolder', required=True)

    identify_parser = subparsers.add_parser('identify', help='Identify GUIDE-seq offtargets')
    identify_parser.add_argument('--aligned', required=True)
    identify_parser.add_argument('--genome', required=True)
    identify_parser.add_argument('--outfolder', required=True)
    identify_parser.add_argument('--target_sequence', required=True)
    identify_parser.add_argument('--description', required=False)
    identify_parser.add_argument('--max_score', required=False, type=int, default=7)
    identify_parser.add_argument('--search_radius', required=False, type=int, default=25)

    filter_parser = subparsers.add_parser('filter', help='Filter identified sites from control sites')
    filter_parser.add_argument('--bedtools', required=True)
    filter_parser.add_argument('--identified', required=True)
    filter_parser.add_argument('--background', required=True)
    filter_parser.add_argument('--outfolder', required=True)

    visualize_parser = subparsers.add_parser('visualize', help='Visualize off-target sites')
    visualize_parser.add_argument('--infile', required=True)
    visualize_parser.add_argument('--outfolder', required=True)
    visualize_parser.add_argument('--title', required=False)

    return parser.parse_args()


def main():
    args = parse_args()

    if args.command == 'all':

        if args.identifyAndFilter:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest)

                # Bootstrap the aligned samfile paths
                g.aligned = {}
                for sample in g.samples:
                    g.aligned[sample] = os.path.join(g.output_folder, 'aligned', sample + '.sam')

                g.identifyOfftargetSites()
                g.filterBackgroundSites()
                g.visualize()

            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_demultiplex:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest)
                g.demultiplexed = {}
                for sample in g.samples:
                    g.demultiplexed[sample] = {}
                    g.demultiplexed[sample]['read1'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.r1.fastq')
                    g.demultiplexed[sample]['read2'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.r2.fastq')
                    g.demultiplexed[sample]['index1'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.i1.fastq')
                    g.demultiplexed[sample]['index2'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.i2.fastq')
                    if not os.path.isfile(g.demultiplexed[sample]['read1']):
                        print ("Can't find ",g.demultiplexed[sample]['read1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['read2']):
                        print ("Can't find ",g.demultiplexed[sample]['read2'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index1']):
                        print ("Can't find ",g.demultiplexed[sample]['index1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index2']):
                        print ("Can't find ",g.demultiplexed[sample]['index2'])
                        exit()

                # Bootstrap the aligned samfile paths
                # g.aligned = {}
                # for sample in g.samples:
                    # g.aligned[sample] = os.path.join(g.output_folder, 'aligned', sample + '.sam')


                g.umitag()
                g.consolidate()
                g.alignReads()
                g.identifyOfftargetSites()
                g.filterBackgroundSites()
                g.visualize()

            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        else:
            g = GuideSeq()
            g.parseManifest(args.manifest)
            g.demultiplex()
            g.umitag()
            g.consolidate()
            g.alignReads()
            g.identifyOfftargetSites()
            g.filterBackgroundSites()
            g.visualize()

    elif args.command == 'demultiplex':
        """
        Run just the demultiplex step given the manifest
        """
        g = GuideSeq()
        g.parseManifestDemultiplex(args.manifest)
        g.demultiplex()

    elif args.command == 'umitag':
        """
        Run just the umitag step
        python guideseq/guideseq.py umitag --read1 test/data/demultiplexed/EMX1.r1.fastq --read2 test/data/demultiplexed/EMX1.r2.fastq --index1 test/data/demultiplexed/EMX1.i1.fastq --index2 test/data/demultiplexed/EMX1.i2.fastq --outfolder test/output/
        """
        g = GuideSeq()
        g.output_folder = args.outfolder
        sample = os.path.basename(args.read1).split('.')[0]
        g.samples = [sample]
        g.demultiplexed = {sample: {}}
        g.demultiplexed[sample]['read1'] = args.read1
        g.demultiplexed[sample]['read2'] = args.read2
        g.demultiplexed[sample]['index1'] = args.index1
        g.demultiplexed[sample]['index2'] = args.index2
        g.umitag()

    elif args.command == 'consolidate':
        """
        Run just the consolidate step
        python guideseq/guideseq.py consolidate --read1 test/data/umitagged/EMX1.r1.umitagged.fastq --read2 test/data/umitagged/EMX1.r2.umitagged.fastq --outfolder test/output/ --min_frequency 0.8 --min_quality 14
        """
        sample = os.path.basename(args.read1).split('.')[0]
        g = GuideSeq()
        g.output_folder = args.outfolder
        g.samples = [sample]
        g.umitagged = {sample: {}}
        g.umitagged[sample]['read1'] = args.read1
        g.umitagged[sample]['read2'] = args.read2

        if 'min_quality' in args:
            min_qual = args.min_quality
        else:
            min_qual = CONSOLIDATE_MIN_QUAL

        if 'min_frequency' in args:
            min_freq = args.min_frequency
        else:
            min_freq = CONSOLIDATE_MIN_FREQ

        g.consolidate(min_freq=min_freq, min_qual=min_qual)

    elif args.command == 'align':
        """
        Run just the alignment step
        python guideseq/guideseq.py align --bwa bwa --read1 test/data/consolidated/EMX1.r1.consolidated.fastq --read2 test/data/consolidated/EMX1.r2.consolidated.fastq --genome /Volumes/Media/hg38/hg38.fa --outfolder test/output/
        """
        sample = os.path.basename(args.read1).split('.')[0]
        g = GuideSeq()
        g.BWA_path = args.bwa
        g.reference_genome = args.genome
        g.output_folder = args.outfolder
        g.samples = [sample]
        g.consolidated = {sample: {}}
        g.consolidated[sample]['read1'] = args.read1
        g.consolidated[sample]['read2'] = args.read2
        g.alignReads()

    elif args.command == 'identify':
        """
        Run just the identify step
        python guideseq/guideseq.py identify --genome /Volumes/Media/hg38/hg38.fa --aligned test/output/aligned/EMX1.sam --outfolder test/output/ --target_sequence GAGTCCGAGCAGAAGAAGAANGG
        """
        if 'description' in args:
            description = args.description
        else:
            description = ''

        if 'max_score' in args:
            max_score = args.max_score
        else:
            max_score = 7

        if 'search_radius' in args:
            search_radius = args.search_radius
        else:
            search_radius = 25

        g = GuideSeq()
        g.output_folder = args.outfolder
        g.reference_genome = args.genome
        sample = os.path.basename(args.aligned).split('.')[0]
        g.samples = {sample: {'description': description, 'target': args.target_sequence}}
        g.aligned = {sample: args.aligned}
        g.max_score = max_score
        g.search_radius = search_radius
        g.identifyOfftargetSites()

    elif args.command == 'filter':
        """
        Run just the filter step

        """
        sample = os.path.basename(args.identified).split('.')[0]
        g = GuideSeq()
        g.output_folder = args.outfolder
        g.bedtools = args.bedtools
        g.samples = {sample: {}, 'control': {}}
        g.identified = {}
        g.identified[sample] = args.identified
        g.identified['control'] = args.background
        g.filterBackgroundSites()

    elif args.command == 'visualize':
        """
        Run just the visualize step
        """
        g = GuideSeq()
        g.output_folder = os.path.dirname(args.outfolder)
        sample = os.path.basename(args.infile).split('.')[0]
        g.samples = {sample: {}}
        g.identified = {}
        g.identified[sample] = args.infile
        g.visualize()


if __name__ == '__main__':
    main()
