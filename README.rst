===============================
GUIDE-Seq
===============================

.. image:: https://img.shields.io/travis/vedtopkar/guideseq.svg
        :target: https://travis-ci.org/vedtopkar/guideseq

.. image:: https://img.shields.io/pypi/v/guideseq.svg
        :target: https://pypi.python.org/pypi/guideseq

.. image:: https://readthedocs.org/projects/guideseq/badge/?version=latest
        :target: https://readthedocs.org/projects/guideseq/?badge=latest
        :alt: Documentation Status


An easy to use bioinformatic pipeline for the GUIDE-seq assay.

* Free software: BSD license
* Documentation: https://guideseq.readthedocs.org.

Features
--------

* TODO


The Pipeline
--------

We start with the 4 FastQ files. The demultiplixer splits these into FastQs for each sample. Consolidation consolidates by molecular index. Aligning is BWA. Then we use the identifyofftarget.py to identify all cleavage sites. Then we bedtools subtract from the control sample to find the actual offtargets (filtering). Finally, we do read sorting and read metrics.
