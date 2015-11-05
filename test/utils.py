import os
import sys
import inspect
import filecmp
from itertools import islice

def checkFolderEquality(folder1, folder2):
    """
    Given two folders, check if there are the same number of files,
    that the names of files are the same, and that the files with the same
    names are the same.
    """

    folder1_files = [x for x in os.listdir(folder1) if not x.startswith('.')]
    folder2_files = [x for x in os.listdir(folder2) if not x.startswith('.')]

    if set(folder1_files) != set(folder2_files):
        print 'Folders do not have the same filenames.'
        return False

    for f in folder1_files:
        if not filecmp.cmp(os.path.join(folder1, f), os.path.join(folder2, f)):
            print '{0} does not match between folders.'.format(f)
            return False

    return True

def head(filepath, n=10):
    with open(filepath) as f:
        for line in islice(f, n):
            print line
