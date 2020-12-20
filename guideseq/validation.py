"""
validation.py
=============

Contains utils for validating the filetype and existence of manifest-defined files/folders

"""

import logging
import os
import sys
from distutils.spawn import find_executable

logger = logging.getLogger('root')


def exists(filepath):
    if not os.path.isfile(filepath):
        logger.error('{0} does not exist'.format(filepath))
        sys.exit()


def checkIfBinary(filepath):
    executable = find_executable(filepath)

    if executable is None:
        logger.error('Executable binary not found at {0}'.format(filepath))
        sys.exit()

    # First check if file exists
    exists(executable)

    # Check if file is a valid binary
    # Adapted from http://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

    if not is_binary_string(open(executable, 'rb').read(1024)):
        logger.error('{0} is not a valid binary'.format(executable))
        sys.exit()


def checkIfFasta(filepath):
    # First check if file exists
    exists(os.path.abspath(filepath))


def checkIfFolder(folderpath):
    # Check if the folder exists
    if not os.path.isdir(os.path.abspath(folderpath)):
        logger.error('{0} is not a valid folder path'.format(folderpath))
        sys.exit()


def checkIfValidUndemultiplexed(undemultiplexed):
    # Check if read1, read2, index1, and index2 exist
    fields = ['forward', 'reverse', 'index1', 'index2']

    if set(fields) != set(undemultiplexed.keys()):
        logger.error('Undemultiplexed field must contain references to "forward", "reverse", "index1", "index2"')
        sys.exit()

    invalid_file = False
    for field in fields:
        if not os.path.isfile(undemultiplexed[field]):
            logger.error('"read1" undemultiplexed field does not reference a valid file')
            invalid_file = True

    if invalid_file:
        sys.exit()


def checkIfValidSamples(samples):
    # Check if control is one of the samples
    if 'control' not in samples:
        logger.error('A control sample must be specified')
        sys.exit()

    if len(samples.keys()) == 0:
        logger.error('No samples defined')
        sys.exit()

    for sample in samples:
        if 'barcode1' not in samples[sample] or 'barcode2' not in samples[sample]:
            logger.error('barcode1 and barcode2 must be specified for {0} sample'.format(sample))
            sys.exit()
        if 'target' not in samples[sample]:
            logger.error('target sequence must be specified for {0} sample'.format(sample))
            sys.exit()


def validateManifest(manifest_data):
    # Check if manifest contains the required fields
    fields = ['bwa', 'bedtools', 'reference_genome', 'output_folder', 'samples', 'undemultiplexed']
    missing_fields = False

    for field in fields:
        if field not in manifest_data.keys():
            logger.error('"{0}" field must be specified in manifest'.format(field))
            missing_fields = True

    if missing_fields:
        sys.exit()

    # Now validate each field
    checkIfBinary(manifest_data['bwa'])
    checkIfBinary(manifest_data['bedtools'])
    checkIfFasta(manifest_data['reference_genome'])
    checkIfValidUndemultiplexed(manifest_data['undemultiplexed'])
    checkIfValidSamples(manifest_data['samples'])
