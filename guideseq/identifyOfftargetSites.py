# IdentifyOffTargetSiteSequences.py
#
# 2015-10-05 Replaced swalign with regex matching.
# 2017-05-31 Replaced nwalign with a explicit search of realignments which uses regex.search.
# 2017-06-03 Output the best offtarget sequences with and/ot without bulges, if any.

from __future__ import print_function

import argparse
import collections
import numpy
import os
import string
import operator
import pyfaidx
import re
import regex
import logging

logger = logging.getLogger('root')


# chromosomePosition defines a class to keep track of the positions.
class chromosomePosition():

	def __init__(self, reference_genome):
		self.chromosome_dict = {}
		self.chromosome_barcode_dict = {}
		self.position_summary = []
		self.index_stack = {}		   # we keep track of the values by index here
		self.genome = pyfaidx.Fasta(reference_genome)

	def addPositionBarcode(self, chromosome, position, strand, barcode, primer, count):
		# Create the chromosome keyValue if it doesn't exist
		if chromosome not in self.chromosome_barcode_dict:
			self.chromosome_barcode_dict[chromosome] = {}
		# Increment the position on that chromosome if it exists, otherwise initialize it with 1
		if position not in self.chromosome_barcode_dict[chromosome]:
			self.chromosome_barcode_dict[chromosome][position] = {}
			self.chromosome_barcode_dict[chromosome][position]['+_total'] = 0
			self.chromosome_barcode_dict[chromosome][position]['+primer1_total'] = 0
			self.chromosome_barcode_dict[chromosome][position]['+primer2_total'] = 0
			self.chromosome_barcode_dict[chromosome][position]['+nomatch_total'] = 0

			self.chromosome_barcode_dict[chromosome][position]['-_total'] = 0
			self.chromosome_barcode_dict[chromosome][position]['-primer1_total'] = 0
			self.chromosome_barcode_dict[chromosome][position]['-primer2_total'] = 0
			self.chromosome_barcode_dict[chromosome][position]['-nomatch_total'] = 0

			self.chromosome_barcode_dict[chromosome][position]['+'] = collections.Counter()
			self.chromosome_barcode_dict[chromosome][position]['+primer1'] = collections.Counter()
			self.chromosome_barcode_dict[chromosome][position]['+primer2'] = collections.Counter()
			self.chromosome_barcode_dict[chromosome][position]['+nomatch'] = collections.Counter()

			self.chromosome_barcode_dict[chromosome][position]['-'] = collections.Counter()
			self.chromosome_barcode_dict[chromosome][position]['-primer1'] = collections.Counter()
			self.chromosome_barcode_dict[chromosome][position]['-primer2'] = collections.Counter()
			self.chromosome_barcode_dict[chromosome][position]['-nomatch'] = collections.Counter()

		self.chromosome_barcode_dict[chromosome][position][strand][barcode] += count
		self.chromosome_barcode_dict[chromosome][position][strand + primer][barcode] += count
		self.chromosome_barcode_dict[chromosome][position][strand + primer + '_total'] += count
		self.chromosome_barcode_dict[chromosome][position][strand + '_total'] += count

	def getSequence(self, genome, chromosome, start, end, strand="+"):
		if strand == "+":
			seq = self.genome[chromosome][int(start):int(end)]
		elif strand == "-":
			seq = self.genome[chromosome][int(start):int(end)].reverse.complement
		return seq

	# Generates a summary of the barcodes by position
	def SummarizeBarcodePositions(self):
		self.barcode_position_summary = [[chromosome, position,
										  len(self.chromosome_barcode_dict[chromosome][position]['+']),
										  len(self.chromosome_barcode_dict[chromosome][position]['-']),
										  self.chromosome_barcode_dict[chromosome][position]['+_total'],
										  self.chromosome_barcode_dict[chromosome][position]['-_total'],
										  len(self.chromosome_barcode_dict[chromosome][position]['+primer1']),
										  len(self.chromosome_barcode_dict[chromosome][position]['+primer2']),
										  len(self.chromosome_barcode_dict[chromosome][position]['-primer1']),
										  len(self.chromosome_barcode_dict[chromosome][position]['-primer2']),
										  ]
										 for chromosome in sorted(self.chromosome_barcode_dict)
										 for position in sorted(self.chromosome_barcode_dict[chromosome])]
		return self.barcode_position_summary

	# Summarizes the chromosome, positions within a 10 bp window
	def SummarizeBarcodeIndex(self, windowsize):
		last_chromosome, last_position, window_index = 0, 0, 0
		index_summary = []
		for chromosome, position, barcode_plus_count, barcode_minus_count, total_plus_count, total_minus_count, plus_primer1_count, plus_primer2_count,\
				minus_primer1_count, minus_primer2_count in self.barcode_position_summary:
			if chromosome != last_chromosome or abs(position - last_position) > 10:
				window_index += 1   # new index
			last_chromosome, last_position = chromosome, position
			if window_index not in self.index_stack:
				self.index_stack[window_index] = []
			self.index_stack[window_index].append([chromosome, int(position),
												   int(barcode_plus_count), int(barcode_minus_count),
												   int(barcode_plus_count) + int(barcode_minus_count),
												   int(total_plus_count), int(total_minus_count),
												   int(total_plus_count) + int(total_minus_count),
												   int(plus_primer1_count), int(plus_primer2_count),
												   int(minus_primer1_count), int(minus_primer2_count)
												   ])
		for index in self.index_stack:
			sorted_list = sorted(self.index_stack[index], key=operator.itemgetter(4))   # sort by barcode_count_total
			chromosome_list, position_list, \
				barcode_plus_count_list, barcode_minus_count_list, barcode_sum_list,\
				total_plus_count_list, total_minus_count_list, total_sum_list, \
				plus_primer1_list, plus_primer2_list, minus_primer1_list, minus_primer2_list\
				= zip(*sorted_list)
			barcode_plus = sum(barcode_plus_count_list)
			barcode_minus = sum(barcode_minus_count_list)
			total_plus = sum(total_plus_count_list)
			total_minus = sum(total_minus_count_list)
			plus_primer1 = sum(plus_primer1_list)
			plus_primer2 = sum(plus_primer2_list)
			minus_primer1 = sum(minus_primer1_list)
			minus_primer2 = sum(minus_primer2_list)
			position_std = numpy.std(position_list)
			min_position = min(position_list)
			max_position = max(position_list)
			barcode_sum = barcode_plus + barcode_minus
			barcode_geometric_mean = (barcode_plus * barcode_minus) ** 0.5
			total_sum = total_plus + total_minus
			total_geometric_mean = (total_plus * total_minus) ** 0.5
			primer1 = plus_primer1 + minus_primer1
			primer2 = plus_primer2 + minus_primer2
			primer_geometric_mean = (primer1 * primer2) ** 0.5
			most_frequent_chromosome = sorted_list[-1][0]
			most_frequent_position = sorted_list[-1][1]
			BED_format_chromosome = "chr" + most_frequent_chromosome
			BED_name = BED_format_chromosome + "_" + str(most_frequent_position) + "_" + str(barcode_sum)
			offtarget_sequence = self.getSequence(self.genome, most_frequent_chromosome, most_frequent_position - windowsize, most_frequent_position + windowsize)

			summary_list = [str(x) for x in [index, most_frequent_chromosome, most_frequent_position, offtarget_sequence,						# pick most frequently occurring chromosome and position
											 BED_format_chromosome, min_position, max_position, BED_name,
											 barcode_plus, barcode_minus, barcode_sum, barcode_geometric_mean,
											 total_plus, total_minus, total_sum, total_geometric_mean,
											 primer1, primer2, primer_geometric_mean, position_std]]

			if (barcode_geometric_mean > 0 or primer_geometric_mean > 0):
				index_summary.append(summary_list)
		return index_summary	# WindowIndex, Chromosome, Position, Plus.mi, Minus.mi,
		# BidirectionalArithmeticMean.mi, BidirectionalGeometricMean.mi,
		# Plus, Minus,
		# BidirectionalArithmeticMean, BidirectionalGeometricMean,


def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
	seq = seq.upper()
	"""
	Given a sequence with ambiguous base characters, returns a regex that matches for
	the explicit (unambiguous) base characters
	"""
	from Bio.Data import IUPACData
	# IUPAC_notation_regex = {'N': '[ATCGN]',
							# 'Y': '[CTY]',
							# 'R': '[AGR]',
							# 'W': '[ATW]',
							# 'S': '[CGS]',
							# 'A': 'A',
							# 'T': 'T',
							# 'C': 'C',
							# 'G': 'G'}
	IUPAC_notation_regex = IUPACData.ambiguous_dna_values
	pattern = ''

	for c in seq:
		pattern += "["+IUPAC_notation_regex[c]+"]"

	if lookahead:
		pattern = '(?b:' + pattern + ')'

	pattern_standard = pattern + '{{s<={0}}}'.format(errors)
	pattern_gap = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
	return pattern_standard, pattern_gap

"""
Realigned TargetSequence and OffTargetSequence when indels are present in the local alignment
"""
def extendedPattern(seq, indels=1, errors=7):
	IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
	realign_pattern = ''
	for c in seq:
		realign_pattern += IUPAC_notation_regex_extended[c]
	return '(?b:' + realign_pattern + ')' + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)


def realignedSequences(targetsite_sequence, chosen_alignment, errors=7):
	match_sequence = chosen_alignment.group()
	substitutions, insertions, deletions = chosen_alignment.fuzzy_counts

	# get the .fuzzy_counts associated to the matching sequence after adjusting for indels, where 0 <= INS, DEL <= 1
	realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))

	if insertions:  # DNA-bulge
		if targetsite_sequence.index('N') > len(targetsite_sequence)/2:  # PAM is on the right end
			targetsite_realignments = [targetsite_sequence[:i + 1] + '-' + targetsite_sequence[i + 1:] for i in range(targetsite_sequence.index('N') + 1)]
		else:
			targetsite_realignments = [targetsite_sequence[:i] + '-' + targetsite_sequence[i:] for i in range(targetsite_sequence.index('N'), len(targetsite_sequence))]
	else:
		targetsite_realignments = [targetsite_sequence]

	realigned_target_sequence, realigned_offtarget_sequence = None, ''  # in case the matching sequence is not founded

	for seq in targetsite_realignments:
		if deletions:  # RNA-bulge
			match_realignments = [match_sequence[:i + 1] + '-' + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
			match_pattern = [match_sequence[:i + 1] + seq[i + 1] + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
		else:
			match_realignments = match_pattern = [match_sequence]

		x = extendedPattern(seq, errors)
		for y_pattern, y_alignment in zip(match_pattern, match_realignments):
			m = regex.search(x, y_pattern, regex.BESTMATCH)
			if m and m.fuzzy_counts == realigned_fuzzy:
				realigned_target_sequence, realigned_offtarget_sequence = seq, y_alignment
	return realigned_target_sequence, realigned_offtarget_sequence


"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match.
"""
def alignSequences(targetsite_sequence, window_sequence, max_score=7):

	window_sequence = window_sequence.upper()
	query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)

	# Try both strands
	alignments_mm, alignments_bulge = list(), list()
	alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
	alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
	alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
	alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

	lowest_distance_score, lowest_mismatch = 100, max_score + 1
	chosen_alignment_b, chosen_alignment_m, chosen_alignment_strand_b, chosen_alignment_strand_m = None, None, '', ''

	for aln_m in alignments_mm:
		strand_m, alignment_type_m, match_m = aln_m
		if match_m != None:
			mismatches, insertions, deletions = match_m.fuzzy_counts
			if mismatches < lowest_mismatch:
				chosen_alignment_m = match_m
				chosen_alignment_strand_m = strand_m
				lowest_mismatch = mismatches

	for aln_b in alignments_bulge:
		strand_b, alignment_type_b, match_b = aln_b
		if match_b != None:
			substitutions, insertions, deletions = match_b.fuzzy_counts
			if insertions or deletions:
				distance_score = substitutions + (insertions + deletions) * 3
				edistance = substitutions + insertions + deletions
				if distance_score < lowest_distance_score and edistance < lowest_mismatch:
					chosen_alignment_b = match_b
					chosen_alignment_strand_b = strand_b
					lowest_distance_score = distance_score

	if chosen_alignment_m:
		offtarget_sequence_no_bulge = chosen_alignment_m.group()
		mismatches = chosen_alignment_m.fuzzy_counts[0]
		start_no_bulge = chosen_alignment_m.start()
		end_no_bulge = chosen_alignment_m.end()
	else:
		offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''

	bulged_offtarget_sequence, score, length, substitutions, insertions, deletions, bulged_start, bulged_end, realigned_target = \
		'', '', '', '', '', '', '', '', 'none'
	if chosen_alignment_b:
		realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_b, max_score)
		if bulged_offtarget_sequence:
			length = len(chosen_alignment_b.group())
			substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
			score = substitutions + (insertions + deletions) * 3
			bulged_start = chosen_alignment_b.start()
			bulged_end = chosen_alignment_b.end()
		else:
			chosen_alignment_strand_b = ''

	return [offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m, start_no_bulge, end_no_bulge,
			bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end, realigned_target]


"""
annotation is in the format:
"""
def analyze(sam_filename, reference_genome, outfile, annotations, windowsize, max_score):

	output_folder = os.path.dirname(outfile)
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	logger.info("Processing SAM file %s", sam_filename)
	file = open(sam_filename, 'rU')
	__, filename_tail = os.path.split(sam_filename)
	chromosome_position = chromosomePosition(reference_genome)
	for line in file:
		fields = line.split('\t')
		if len(fields) >= 10:
			# These are strings--need to be cast as ints for comparisons.
			full_read_name, sam_flag, chromosome, position, mapq, cigar, name_of_mate, position_of_mate, template_length, read_sequence, read_quality = fields[:11]
			if int(mapq) >= 50 and int(sam_flag) & 128 and not int(sam_flag) & 2048:
				# Second read in pair
				barcode, count = parseReadName(full_read_name)
				primer = assignPrimerstoReads(read_sequence, sam_flag)
				if int(template_length) < 0:  # Reverse read
					read_position = int(position_of_mate) + abs(int(template_length)) - 1
					strand = "-"
					chromosome_position.addPositionBarcode(chromosome, read_position, strand, barcode, primer, count)
				elif int(template_length) > 0:  # Forward read
					read_position = int(position)
					strand = "+"
					chromosome_position.addPositionBarcode(chromosome, read_position, strand, barcode, primer, count)

	# Generate barcode position summary
	stacked_summary = chromosome_position.SummarizeBarcodePositions()

	with open(outfile, 'w') as f:
		# Write header
		print('#BED_Chromosome', 'BED_Min.Position', 'BED_Max.Position', 'BED_Name', 'Filename',
			  'WindowIndex', 'WindowChromosome', 'Position', 'WindowSequence',
			  '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi', '+.total', '-.total', 'total.sum', 'total.geometric_mean',
			  'primer1.mi', 'primer2.mi', 'primer.geometric_mean', 'position.stdev',
			  'BED_Site_Name', 'BED_Score', 'BED_Site_Chromosome',
			  'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',  # 24:25
			  'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',  # 26:28
			  'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',  # 29:31
			  'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',  # 32:34
			  'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',  # 35:37
			  'Cell', 'Targetsite', 'TargetSequence', 'RealignedTargetSequence', sep='\t', file=f)  # 38:41

		# Output summary of each window
		summary = chromosome_position.SummarizeBarcodeIndex(windowsize)
		target_sequence = annotations["Sequence"]
		annotation = [annotations['Description'],
					  annotations['Targetsite'],
					  annotations['Sequence']]
		output_dict = {}

		for row in summary:
			window_sequence, window_chromosome, window_start, window_end, BED_name = row[3:8]

			non_bulged_target_start_absolute, bulged_target_start_absolute = '', ''
			if target_sequence:
				offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
				bulged_offtarget_sequence, length, distance, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end, \
				realigned_target_sequence = alignSequences(target_sequence, window_sequence, max_score)

				BED_score = 1
				BED_chromosome = window_chromosome
				if chosen_alignment_strand_m == "+":
					non_bulged_target_start_absolute = start_no_bulge + int(row[2]) - windowsize
					non_bulged_target_end_absolute = end_no_bulge + int(row[2]) - windowsize
				elif chosen_alignment_strand_m == "-":
					non_bulged_target_start_absolute = int(row[2]) + windowsize - end_no_bulge
					non_bulged_target_end_absolute = int(row[2]) + windowsize - start_no_bulge
				else:
					non_bulged_target_start_absolute, non_bulged_target_end_absolute = [""] * 2

				if chosen_alignment_strand_b == "+":
					bulged_target_start_absolute = bulged_start + int(row[2]) - windowsize
					bulged_target_end_absolute = bulged_end + int(row[2]) - windowsize
				elif chosen_alignment_strand_b == "-":
					bulged_target_start_absolute = int(row[2]) + windowsize - bulged_end
					bulged_target_end_absolute = int(row[2]) + windowsize - bulged_start
				else:
					bulged_target_start_absolute, bulged_target_end_absolute = [""] * 2

				if not (chosen_alignment_strand_m or chosen_alignment_strand_b):
					BED_chromosome, BED_score, BED_name = [""] * 3
				# print ("BED_name",BED_name)
				# print ("annotation",annotation)
				output_row = row[4:8] + [filename_tail] + row[0:4] + row[8:] + \
							 [str(BED_name), BED_score, BED_chromosome,
											  offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m,
											  non_bulged_target_start_absolute, non_bulged_target_end_absolute,
											  bulged_offtarget_sequence, length, distance, substitutions, insertions, deletions,
											  chosen_alignment_strand_b, bulged_target_start_absolute, bulged_target_end_absolute] + \
							 [str(x) for x in annotation] + [realigned_target_sequence]
			else:
				output_row = [str(x) for x in row[4:8] + [filename_tail] + row[0:4] + row[8:] + [""] * 17 + annotation + ['none']]

			if non_bulged_target_start_absolute != '' or bulged_target_start_absolute != '':
				# print ("non_bulged_target_start_absolute",non_bulged_target_start_absolute)
				# print ("bulged_target_start_absolute",bulged_target_start_absolute)
				# print ("non_bulged_target_end_absolute",non_bulged_target_end_absolute)
				# print ("bulged_target_end_absolute",bulged_target_end_absolute)
				output_row_key = '{0}_{1}_{2}'.format(window_chromosome, py2min([non_bulged_target_start_absolute, bulged_target_start_absolute]), py2max([non_bulged_target_end_absolute, bulged_target_end_absolute]))
			else:
				output_row_key = '{0}_{1}_{2}'.format(window_chromosome, window_start, window_end)

			if output_row_key in output_dict.keys():
				read_count_total = int(output_row[11]) + int(output_dict[output_row_key][11])
				output_dict[output_row_key][11] = str(read_count_total)
			else:
				output_dict[output_row_key] = output_row

		for key in sorted(output_dict.keys()):
			print(*output_dict[key], sep='\t', file=f)

def py2min(myList):
	out = [i  for i in myList if i != "" ]
	return min(out)
	
def py2max(myList):
	out = [i  for i in myList if i != "" ]
	return max(out)

def assignPrimerstoReads(read_sequence, sam_flag):
	# Get 20-nucleotide sequence from beginning or end of sequence depending on orientation
	if int(sam_flag) & 16:
		readstart = reverseComplement(read_sequence[-20:])
	else:
		readstart = read_sequence[:20]
	if readstart == "TTGAGTTGTCATATGTTAAT":
		return "primer1"
	elif readstart == "ACATATGACAACTCAATTAA":
		return "primer2"
	else:
		return "nomatch"


def loadFileIntoArray(filename):
	with open(filename, 'rU') as f:
		keys = f.readline().rstrip('\r\n').split('\t')[1:]
		data = collections.defaultdict(dict)
		for line in f:
			filename, rest = processLine(line)
			line_to_dict = dict(zip(keys, rest))
			data[filename] = line_to_dict
	return data


def parseReadName(read_name):
	m = re.search(r'([ACGTN]{8}_[ACGTN]{6}_[ACGTN]{6})_([0-9]*)', read_name)
	if m:
		molecular_index, count = m.group(1), m.group(2)
		return molecular_index, int(count)
	else:
		# print read_name
		return None, None


def processLine(line):
	fields = line.rstrip('\r\n').split('\t')
	filename = fields[0]
	rest = fields[1:]
	return filename, rest

import sys
def reverseComplement(sequence):
	if sys.version_info[0] < 3:
		tab = string.maketrans("ACGTacgt", "TGCATGCA")
	else:
		tab = str.maketrans("ACGTacgt", "TGCATGCA")
	return sequence.translate(tab)[::-1]


def main():
	parser = argparse.ArgumentParser(description='Identify off-target candidates from Illumina short read sequencing data.')
	parser.add_argument('--ref', help='Reference Genome Fasta', required=True)
	parser.add_argument('--samfile', help='SAM file', nargs='*')
	parser.add_argument('--outfile', help='File to output identified sites to.', required=True)
	parser.add_argument('--window', help='Window around breakpoint to search for off-target', type=int, default=25)
	parser.add_argument('--max_score', help='Score threshold', type=int, default=7)
	# parser.add_argument('--demo')
	parser.add_argument('--target', default='')

	args = parser.parse_args()

	annotations = {'Description': 'test description', 'Targetsite': 'dummy targetsite', 'Sequence': args.target}
	analyze(args.samfile[0], args.ref, args.outfile, annotations, args.window, args.max_score)


if __name__ == "__main__":
	main()
