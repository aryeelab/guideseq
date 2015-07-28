#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_guideseq
----------------------------------

Tests for `guideseq` module.
"""

import unittest

from guideseq import guideseq

class DemultiplexTestCase(unittest.TestCase):

    def setUp(self):
        # do demultiplexing
        pass

    def testIfDemultiplexed(self):
        # See if equal to known demultiplexed output
        pass

    def tearDown(self):
        # Delete temp output
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