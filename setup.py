#! /usr/bin/python
# coding=utf-8

from setuptools import setup, find_packages

if __name__ == '__main__':
    RELEASE = "1.2.0.dev"
    with open("README.rst") as summary:
        LONG_DESCRIPTION = summary.read()
    CLASSIFIERS = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
    setup(name='RotamerConvolveMD',
          version=RELEASE,
          description='Analysis of spin label distances over structural ensembles',
          author='Philip W. Fowler',
          author_email='philip.fowler@ndm.ox.ac.uk',
          url='https://github.com/MDAnalysis/RotamerConvolveMD',
          requires=['numpy (>=1.6)', 'MDAnalysis (>=0.11.0)'],
          install_requires=['numpy>=1.6.0', 'MDAnalysis>=0.11.0'],
          provides=['rotcon'],
          license='GPL 2',
          packages=find_packages(exclude=['scripts', 'rotcon/data']),
          package_data={'rotcon': ['data/*.pdb', 'data/*.dcd', 'data/*.dat']},
          scripts=['scripts/convolve-mtss-rotamers.py', 'scripts/convolve-mtss-rotamers_pre.py'],
          classifiers=CLASSIFIERS,
          long_description=LONG_DESCRIPTION,
          zip_safe=True,
          )
