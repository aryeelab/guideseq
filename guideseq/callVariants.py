from __future__ import print_function

import subprocess
import sys
import os
import argparse
import regex
import re
import HTSeq
import pyfaidx
from findCleavageSites import get_sequence, regexFromSequence, alignSequences, reverseComplement, extendedPattern, realignedSequences


"""
Run samtools:mpileup and get all identified variants in the window sequences
"""
def snpCall(matched_file, reference, bam_file, out, search_radius):
    basename = os.path.basename(out)
    output_folder = os.path.dirname(out)

    # open matched file
    regions = list()
    with open(matched_file, 'rU') as f:
        f.readline()
        for line in f:
            site = line.strip().split('\t')
            regions.append([site[0], int(site[1]) - search_radius, int(site[2]) + search_radius, '*', bam_file, '_'.join([site[34], site[3]])])

    print('Running samtools:mpileup for %s' % basename, file=sys.stderr)
    out_vcf = os.path.join(output_folder, basename + '_mpileup_output')
    if os.path.exists(out_vcf):
        subprocess.check_call('rm -r %s' % out_vcf, shell=True, env=os.environ.copy())
    os.makedirs(out_vcf)
    process_mpileup = open(os.path.join(out_vcf, 'logFile_mpileup'), 'w')

    for item in regions:
        chromosome, start, end, strand, bam_file, region_basename = item
        region = '%s%s%s%s%s' % (chromosome, ":", int(start), "-", int(end))
        output = os.path.join(out_vcf, region_basename + '.vcf')

        cl_vcf = 'samtools mpileup -v --region %s --fasta-ref %s %s > %s' % (region, reference, bam_file, output)
        subprocess.check_call(cl_vcf, shell=True, env=os.environ.copy(), stderr=process_mpileup, stdout=process_mpileup)
    process_mpileup.close()

    print('Collecting variants for %s' % basename, file=sys.stderr)
    out_bcf = os.path.join(output_folder, basename + '_output_bcftools')
    if os.path.exists(out_bcf):
        subprocess.check_call('rm -r %s' % out_bcf, shell=True, env=os.environ.copy())
    os.makedirs(out_bcf)
    process_bcftools = open(os.path.join(out_bcf, 'logFile_bcftools'), 'w')

    vcf_files = [f for f in os.listdir(out_vcf) if os.path.isfile(os.path.join(out_vcf, f))]
    for arch in vcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf'):
            name = arch[:-4]
            output = os.path.join(out_bcf, name + '_BCFcall.vcf')

            cl_bcf = 'bcftools call -v -c %s > %s' % (os.path.join(out_vcf, arch), output)
            subprocess.check_call(cl_bcf, shell=True, env=os.environ.copy(), stderr=process_bcftools, stdout=process_bcftools)
    process_bcftools.close()

    print('Collecting significant variant calls for %s' % basename, file=sys.stderr)
    out_svc = os.path.join(output_folder, basename + '_output_svc')
    if os.path.exists(out_svc):
        subprocess.check_call('rm -r %s' % out_svc, shell=True, env=os.environ.copy())
    os.makedirs(out_svc)
    process_svc = open(os.path.join(out_svc, 'logFile_svc'), 'w')

    bcf_files = [f for f in os.listdir(out_bcf) if os.path.isfile(os.path.join(out_bcf, f))]
    for arch in bcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf'):
            name = arch[:-12]
            output = os.path.join(out_svc, name + '_SIGNFcall.txt')

            cl_sed = "sed -n '/##/!p' %s | awk 'FNR>1' > %s" % (os.path.join(out_bcf, arch), output)
            subprocess.check_call(cl_sed, shell=True, env=os.environ.copy(), stderr=process_svc, stdout=process_svc)
    process_svc.close()

    print('Consolidating all the significant variant calls for %s' % basename, file=sys.stderr)
    header = ['targetsite', 'site_name', 'chromosome', 'one_based_position', 'reference', 'variant', 'quality', 'genotype', 'depth', 'PL']
    variants = list()

    svc_files = [f for f in os.listdir(out_svc) if os.path.isfile(os.path.join(out_svc, f))]
    for arch in svc_files:
        if not arch.startswith('.') and arch.endswith('.txt'):
            tag = arch[:-14]
            f = open(os.path.join(out_svc, arch), 'r')
            reads = f.readlines()
            f.close()

            for line in reads:
                item = line.split()
                if 'INDEL' in item[7]:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[3][3:]] + ['_'.join(item[9][4:].split(','))])
                else:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[0][3:]] + ['_'.join(item[9][4:].split(','))])

    out_file = open(out + '_mpileupCall.txt', 'w')
    print(*header, sep='\t', file=out_file)
    for item in variants:
        print(*item, sep='\t', file=out_file)
    out_file.close()

    print('Cleaning up directive for %s' % basename, file=sys.stderr)
    subprocess.check_call('rm -r %s' % out_vcf, shell=True, env=os.environ.copy())
    subprocess.check_call('rm -r %s' % out_bcf, shell=True, env=os.environ.copy())
    subprocess.check_call('rm -r %s' % out_svc, shell=True, env=os.environ.copy())

    print('Done running samtools:mpileup for %s' % basename, file=sys.stderr)
    return variants


"""
Obtain variant off-target sequences
"""
def realignVariantBulge(bulge_sequence, window_sequence_variant, bulge_strand):
    bseq = bulge_sequence.replace('-', '')
    if bulge_strand == '+':
        m_bulge = re.search(bseq, window_sequence_variant, re.I)
    else:
        m_bulge = re.search(bseq, reverseComplement(window_sequence_variant), re.I)
    variant_bseq = m_bulge.group()
    variant_bseq = variant_bseq[:bulge_sequence.find('-')] + '-' + variant_bseq[bulge_sequence.find('-'):]
    return variant_bseq


def SNPreader(snp_file):
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode='O')

    for snp in snp_file:
        basename, snpID, chromosome, one_based_position, reference, variant, quality, genotype, depth, PL = snp
        position = int(one_based_position) - 1
        key = '_'.join([basename, chromosome])
        ga[HTSeq.GenomicInterval(chromosome, position, position + 1, ".")] = [position, reference, variant, genotype, quality, key]
    return ga


def arrayOffTargets(matched_file, search_radius):
    offtargets_dict = {}
    gi_dict = {}

    with open(matched_file, 'r') as g:
        g.readline()
        for line in g:
            site = line.strip().split('\t')

            Chromosome = site[0]
            start = int(site[1]) - search_radius
            end = int(site[2]) + search_radius
            Name = site[3]

            offtargets_dict[Name] = site

            gi_dict[Name] = HTSeq.GenomicInterval(Chromosome, start, end, ".")
    return offtargets_dict, gi_dict


def snpAdjustment(matched_file, snp_file, out, mismatch_threshold, search_radius):
    output_file = open(out + '_Variants.txt', 'w')
    print('Chromosome', 'Start', 'End', 'Name', 'ReadCount',
          'Variant_WindowSequence',
          'Variant_Site_SubstitutionsOnly.Sequence', 'Variant_Site_SubstitutionsOnly.NumSubstitutions',
          'Variant_Site_SubstitutionsOnly.Strand',
          'Variant_Site_GapsAllowed.Sequence', 'Variant_Site_GapsAllowed.Length', 'Variant_Site_GapsAllowed.Score',
          'Variant_Site_GapsAllowed.Substitutions', 'Variant_Site_GapsAllowed.Insertions', 'Variant_Site_GapsAllowed.Deletions',
          'Variant_Site_GapsAllowed.Strand',
          'Cell', 'Targetsite', 'TargetSequence', 'Variant_RealignedTargetSequence',
          'Reference', 'Variant', 'Genotype', 'Quality',
          sep='\t', file=output_file)
    output_file.close()

    basename = os.path.basename(out)
    offtargets, gi_offtargets = arrayOffTargets(matched_file, search_radius)
    ga_snp = SNPreader(snp_file)

    for name in offtargets:
        variant_flag = False
        site = offtargets[name]
        gi = gi_offtargets[name]

        chromosome = site[0]
        window_sequence = site[6]
        window_sequence = window_sequence.upper()
        cell, targetsite, TargetSequence = site[33:36]
        output01 = site[0:4] + [site[9]]
        output03 = [cell, targetsite, TargetSequence]
        ots_nb, ots_bu = site[19], site[24]

        #  obtain variant window sequence
        wkey = '_'.join([basename, chromosome])
        insert_start, insert_end, insert_var, snp_data = list(), list(), list(), {}

        for i, v in ga_snp[gi].steps():
            if v:
                position, reference, variant, genotype, quality, key = v
                if key == wkey:
                    variant = variant.split(',')[0]
                    for n, pos in enumerate(range(gi.start, gi.end)):
                        if pos == int(position):
                            insert_var.append(variant.lower())
                            insert_start.append(n)
                            end_pos = n + len(reference)
                            insert_end.append(end_pos)
                            snp_data[str(position)] = [position, reference, variant, genotype, quality]

        tri = 0
        window_sequence_variant = ''
        for i in range(len(insert_var)):
            variant = insert_var[i]
            pos = insert_start[i]
            window_sequence_variant += window_sequence[tri:pos] + variant.lower()
            tri = insert_end[i]
        window_sequence_variant += window_sequence[tri:]

        #  variant off-target sequences: only proceed if there is a variant in the window sequence
        window_sequence_var = window_sequence_variant.upper()
        if window_sequence_var != window_sequence:
            offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, \
            chosen_alignment_strand_b, bulged_start, bulged_end, realigned_target = \
                alignSequences(TargetSequence, window_sequence_var, max_score=mismatch_threshold)

            variant_ots_no_bulge, variant_ots_bulge = '', ''

            #  get variant sequence if the off-target sequences changed by considering the variant window
            if ots_nb != offtarget_sequence_no_bulge:
                variant_flag = True
                if chosen_alignment_strand_m == '+':
                    m_no_bulge = re.search(offtarget_sequence_no_bulge, window_sequence_variant, re.I)
                else:
                    m_no_bulge = re.search(offtarget_sequence_no_bulge, reverseComplement(window_sequence_variant), re.I)
                variant_ots_no_bulge = m_no_bulge.group()

            if ots_bu != bulged_offtarget_sequence:
                variant_flag = True
                variant_ots_bulge = realignVariantBulge(bulged_offtarget_sequence, window_sequence_variant, chosen_alignment_strand_b)

            # collect and write variant data if we have variant off-target sequence(s)
            if variant_flag:
                total_genotype, total_reference, total_variant, total_quality = '', '', '', ''
                for pos in snp_data:
                    position, reference, variant, genotype, quality = snp_data[pos]
                    if total_genotype != '':
                        total_genotype += ''.join([':', genotype])
                        total_reference += ''.join([':', reference])
                        total_variant += ''.join([':', variant])
                        total_quality += ''.join([':', quality])
                    else:
                        total_genotype += ''.join([genotype])
                        total_reference += ''.join([reference])
                        total_variant += ''.join([variant])
                        total_quality += ''.join([quality])

                output02 = [variant_ots_no_bulge, mismatches, chosen_alignment_strand_m,
                            variant_ots_bulge, length, substitutions, insertions, deletions, chosen_alignment_strand_b]
                output04 = [total_reference, total_variant, total_genotype, total_quality]
                output_line = output01 + [window_sequence_variant] + output02 + output03 + [realigned_target] + output04

                with open(out + '_Variants.txt', 'a') as output_file:
                    print(*output_line, sep='\t', file=output_file)

def sortSAM(sam_file):
    sam_folder = os.path.dirname(sam_file)
    sam_name = os.path.basename(sam_file)
    bam_file = os.path.join(sam_folder, sam_name[:-4] + '.bam')

    samtools_sam_to_bam_command = 'samtools sort {0}'.format(bam_file, sam_file)
    samtools_index_command = 'samtools index {0}'.format(bam_file)

    # Convert SAM to BAM file
    subprocess.check_call(samtools_sam_to_bam_command, shell=True)

    # Index BAM file
    subprocess.check_call(samtools_index_command, shell=True)


"""
Main function
"""
def getVariants(matched_file, ref, sam_file, out, search_radius, mismatch_threshold):
    basename = os.path.basename(out)
    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    sortSAM(sam_file)
    bam_file = sam_file[:-4] + '.bam'

    snp_file = snpCall(matched_file, ref, bam_file, out, search_radius)

    print('Obtaining Variant Off-Target Sequences for %s' % basename, file=sys.stderr)
    snpAdjustment(matched_file, snp_file, out, mismatch_threshold, search_radius)


def main():
    parser = argparse.ArgumentParser(description='Implement samtools:mpileup to identify genomic variants and adjust the off-target sequence when required.')
    parser.add_argument('--matched_file', help="full_path_to/matched file in 'identified' folder", required=True)
    parser.add_argument('--ref', help="Reference Genome Fasta", required=True)
    parser.add_argument('--sam', help="SAM file", required=True)
    parser.add_argument('--search_radius', help="Search radius around the position window", default=20, type=int)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=7, type=int)
    parser.add_argument('--out', help="Output file basename, with full path", required=True)
    args = parser.parse_args()

    getVariants(args.matched_file, args.ref, args.sam[0], args.out, args.search_radius, args.mismatch_threshold)

if __name__ == "__main__":
    main()
