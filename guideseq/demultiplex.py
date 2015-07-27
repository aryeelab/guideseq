"""

demultiplex.py

Authored by Martin Aryee, modified by Ved Topkar

"""

from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time

__author__ = 'Martin Aryee'

def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')

    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]


def get_sample_id(i1, i2, sample_names):
    seq1 = i1[1]
    seq2 = i2[1]
    sample_barcode = seq1[1:8] + seq2[1:8]
    if sample_barcode in sample_names:
        return sample_names[sample_barcode]
    else:
        return sample_barcode


def demultiplex(read1, read2, index1, index2, out_dir, sample_barcodes, min_reads=10000):

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    sample_names = dict(zip(sample_barcodes.values(), sample_barcodes.keys()))

    outfiles_r1 = {}
    outfiles_r2 = {}
    outfiles_i1 = {}
    outfiles_i2 = {}

    total_count = 0
    count = {}
    buffer_r1 = {}
    buffer_r2 = {}
    buffer_i1 = {}
    buffer_i2 = {}

    start = time.time()
    for r1,r2,i1,i2 in itertools.izip(fq(read1), fq(read2), fq(index1), fq(index2)):
        total_count += 1
        if total_count % 1000000 == 0:
            print ("Processed %d reads in %.1f minutes." % (total_count, (time.time()-start)/60))
        sample_id = get_sample_id(i1, i2, sample_names)

        # Increment read count and create output buffers if this is a new sample barcode
        if not sample_id in count:
            count[sample_id] = 0
            buffer_r1[sample_id] = []
            buffer_r2[sample_id] = []
            buffer_i1[sample_id] = []
            buffer_i2[sample_id] = []
        count[sample_id] += 1

        # Write the reads to file or to sample-specific output buffers if we haven't
        # reached min_reads
        # The buffers are used to avoid writing fastqs for sample barcodes with very few reads.
        if count[sample_id] < min_reads:
            buffer_r1[sample_id].append(r1)
            buffer_r2[sample_id].append(r2)
            buffer_i1[sample_id].append(i1)
            buffer_i2[sample_id].append(i2)
        elif count[sample_id] == min_reads:
            outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'w')
            outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'w')
            outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % sample_id), 'w')
            outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % sample_id), 'w')
            # Spill the buffers to sample-specific fastqs
            for record in buffer_r1[sample_id] + r1:
                outfiles_r1[sample_id].write(''.join(record))
            for record in buffer_r2[sample_id] + r2:
                outfiles_r2[sample_id].write(''.join(record))
            for record in buffer_i1[sample_id] + i1:
                outfiles_i1[sample_id].write(''.join(record))
            for record in buffer_i2[sample_id] + i2:
                outfiles_i2[sample_id].write(''.join(record))
            del buffer_r1[sample_id]
            del buffer_r2[sample_id]
            del buffer_i1[sample_id]
            del buffer_i2[sample_id]
        else:
            for line in r1:
                print (line, file=outfiles_r1[sample_id], end="")
            for line in r2:
                print (line, file=outfiles_r2[sample_id], end="")
            for line in i1:
                print (line, file=outfiles_i1[sample_id], end="")
            for line in i2:
                print (line, file=outfiles_i2[sample_id], end="")

    # Write remaining buffered reads to a single fastq.
    # (These reads correspond to barcodes that were seen less than min_reads times)
    undetermined_r1 = open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'w')
    undetermined_r2 = open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'w')
    undetermined_i1 = open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'w')
    undetermined_i2 = open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'w')
    for sample_id in buffer_r1.keys():
        for record in buffer_r1[sample_id]:
            undetermined_r1.write(''.join(record))
        for record in buffer_r2[sample_id]:
            undetermined_r2.write(''.join(record))
        for record in buffer_i1[sample_id]:
            undetermined_i1.write(''.join(record))
        for record in buffer_i2[sample_id]:
            undetermined_i2.write(''.join(record))

    # Close files
    for sample_id in outfiles_r1:
        outfiles_r1[sample_id].close()
        outfiles_r2[sample_id].close()
        outfiles_i1[sample_id].close()
        outfiles_i2[sample_id].close()
    undetermined_r1.close()
    undetermined_r2.close()
    undetermined_i1.close()
    undetermined_i2.close()

    num_fastqs = len([v for k,v in count.iteritems() if v >= min_reads])
    print('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.' % (num_fastqs, len(count), min_reads))

    print (outfiles_r1, outfiles_r2, outfiles_i1, outfiles_i2)
    return (outfiles_r1, outfiles_r2, outfiles_i1, outfiles_i2)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--read1', required=True)
    parser.add_argument('--read2', required=True)
    parser.add_argument('--index1', required=True)
    parser.add_argument('--index2', required=True)
    parser.add_argument('--min_reads', type=int, default=10000)
    parser.add_argument('--sample_barcodes')
    parser.add_argument('--out_dir', default='.')

    args = vars(parser.parse_args())

    read1 = args['read1']
    read2 = args['read2']
    index1 = args['index1']
    index2 = args['index2]']

    out_dir = args['out_dir']

    demultiplex(read1, read2, index1, index2, out_dir, min_reads=10000)


if __name__ == '__main__':
    main()
