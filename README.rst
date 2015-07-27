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

* Python (2.6, 2.7, or PyPy)
* `bwa <http://bio-bwa.sourceforge.net/>`_ alignment tool
* `bedtools <http://bedtools.readthedocs.org/en/latest/>`_.
* Reference genome .fasta file (we recommend `hg19 <http://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=431825753_a0WJjTe0PI8wUUlzy80AAMLzPJg4&clade=mammal&org=Human&db=hg19>`_)

Usage
=======

Using this tool is simple, just create a ``.yaml`` manifest file referencing the dependencies and sample ``.fastq.gz`` file paths. Below is an example::

    reference_genome: /Volumes/Media/hg38/hg38.fa
    output_folder: ../tests/output

    bwa: /Users/VedTopkar/code/bwa/bwa
    bedtools: bedtools

    undemultiplexed:
        forward: ../tests/data/undemux.r1.fastq.gz
        reverse: ../tests/data/undemux.r2.fastq.gz
        index1: ../tests/data/undemux.i1.fastq.gz
        index2: ../tests/data/undemux.i2.fastq.gz

    sample_barcodes:
        control: AGGCATGAGATCGC
        EMX1: GACTCCTGCGATAT

Absolute paths are recommended. Be sure to point the ``bwa`` and ``bedtools`` paths directly to their respective executables.

Once you have a manifest file created, you can simply execute ``python guideseq.py -m PATH/TO/MANIFEST.YAML`` to run the entire pipeline.

You cannot run steps of the pipeline individually, though this functionality is planned for future releases.