import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('sim-reads requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = ['suds', 'biopython']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import sim_reads as distmeta

setup(
    name='sim-reads',
    version=distmeta.__version__,
    description='Generate simulated reads.',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['sim_reads'],
    install_requires=requires,
    entry_points = {
        'console_scripts': [
            'make_fastq = sim_reads.make_fastq:main'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
