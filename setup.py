#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

requirements = [line.strip() for line in open('requirements.txt')]

setup(
    name='guideseq',
    version='0.9.0',
    description="An easy to use bioinformatic pipeline for the GUIDE-seq assay.",
    author="Shengdar Q Tsai, Martin Aryee, Ved V Topkar",
    author_email='STSAI4@mgh.harvard.edu, Aryee.Martin@mgh.harvard.edu, vedtopkar@gmail.com',
    url='https://github.com/vedtopkar/guideseq',
    packages=[
        'guideseq',
    ],
    package_dir={'guideseq':
                 'guideseq'},
    install_requires=requirements,
    license="AGPL",
    keywords='guideseq',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: Unix',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7'
    ]
)
