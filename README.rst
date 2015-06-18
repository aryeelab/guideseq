===============================
GUIDE-Seq
===============================

.. image:: https://img.shields.io/travis/vedtopkar/guideseq.svg
        :target: https://travis-ci.org/vedtopkar/guideseq

.. image:: https://coveralls.io/repos/vedtopkar/guideseq/badge.svg?branch=master 
        :target: https://coveralls.io/r/vedtopkar/guideseq?branch=master

.. image:: https://img.shields.io/pypi/v/guideseq.svg
        :target: https://pypi.python.org/pypi/guideseq

.. image:: https://readthedocs.org/projects/guideseq/badge/?version=latest
        :target: http://guideseq.readthedocs.org/en/latest/
        :alt: Documentation Status


GUIDE-Seq An easy to use bioinformatic pipeline for the GUIDE-seq assay.

* Free software: License not yet determined.
* Documentation: https://guideseq.readthedocs.org.



Dependencies
=======

* Python (2.6, 2.7, PyPy)
* `bwa <http://bio-bwa.sourceforge.net/>`_ alignment tool
* `bedtools <http://bedtools.readthedocs.org/en/latest/>`_.
* Reference genome .fasta file (we recommend HG19)

Usage
=======

Using this tool is simple, just create a 


What is GUIDE-Seq?
========

Genome-wide, Unbiased Identification of DSBs Enabled by Sequencing (GUIDE-seq) is a method for ...


The Pipeline
========

.. Demultiplexing
.. --------

.. The demultiplexing step splits the four FASTQ input files into FASTQ files for each sample (as specified in the manifest.yaml file).


.. Consolidation
.. --------

.. The consolidation step consolidates reads by their molecular indices.


.. Alignment
.. --------

.. In this step, the reads are aliged to HG19 (as provided) to map the detected DSBs. This step requires the BWA alignment program.


.. Identify Offtargets

.. Filtering
.. ---------


We start with the 4 FastQ files. The demultiplixer splits these into FastQs for each sample. Consolidation consolidates by molecular index. Aligning is BWA. Then we use the identifyofftarget.py to identify all cleavage sites. Then we bedtools subtract from the control sample to find the actual offtargets (filtering). Finally, we do read sorting and read metrics.