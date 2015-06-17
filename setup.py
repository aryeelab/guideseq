#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    'PyYAML'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='guideseq',
    version='0.1.0',
    description="An easy to use bioinformatic pipeline for the GUIDE-seq assay.",
    long_description=readme + '\n\n' + history,
    author="Shengdar Q Tsai, Martin Aryee, Ved Topkar",
    author_email='STSAI4@mgh.harvard.edu, Aryee.Martin@mgh.harvard.edu, vedtopkar@college.harvard.edu',
    url='https://github.com/vedtopkar/guideseq',
    packages=[
        'guideseq',
    ],
    package_dir={'guideseq':
                 'guideseq'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords='guideseq',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7'
    ],
    test_suite='tests',
    tests_require=test_requirements
)
