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
        file1 = os.path.join(folder1, f)
        file2 = os.path.join(folder2, f)

        if f.split('.')[-1] == 'sam':
            with open(file1, 'r') as a, open(file2, 'r') as b:
                for line1, line2 in zip(a,b):
                    if line1.startswith('@'):
                        continue
                    elif line1 != line2:
                        return False
        else:
            if not filecmp.cmp(file1, file2):
                print '{0} does not match between folders.'.format(f)
                return False

    return True


def head(filepath, n=10):
    with open(filepath) as f:
        for line in islice(f, n):
            print line
