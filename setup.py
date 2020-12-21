#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import guideseq

import os
if os.path.isfile("README.MD"):
	with open("README.MD", "r") as fh:
		long_description = fh.read()
else:
	long_description="guide-seq"


setup(
	name='guide_seq',
	version=str(guideseq.__version__),
	description="An easy to use bioinformatic pipeline for the GUIDE-seq assay.",
	author="Shengdar Q Tsai, Martin Aryee, Ved V Topkar",
	author_email='STSAI4@mgh.harvard.edu, Aryee.Martin@mgh.harvard.edu, vedtopkar@gmail.com',
	url='https://github.com/tsailabSJ/guideseq',
	# packages=find_packages(),
	packages=[
		'guideseq',
		'umi',
	],
	package_dir={'guideseq':
				 'guideseq','umi':'guideseq/umi'},
	
	scripts=['guideseq/guideseq.py','guideseq/alignReads.py','guideseq/visualization.py',
		'guideseq/filterBackgroundSites.py','guideseq/identifyOfftargetSites.py','guideseq/log.py',
		'guideseq/validation.py'],
	package_data={'test': ["test/*"]},
	license="AGPL",
	include_package_data=True,
	long_description=long_description,
	long_description_content_type='text/markdown',
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
		'Programming Language :: Python :: 3'
	]
)
