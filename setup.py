#! /usr/bin/python
# coding=utf-8

from setuptools import setup, find_packages

if __name__ == '__main__':
    RELEASE = "1.3.0"
    with open("README.rst") as summary:
        LONG_DESCRIPTION = summary.read()


    setup(name='RotamerConvolveMD',
          version=RELEASE,
          description='Analysis of spin label distances over structural ensembles',
          author='Philip W. Fowler',
          author_email='philip.fowler@ndm.ox.ac.uk',
          url='https://github.com/MDAnalysis/RotamerConvolveMD',
          requires=['numpy (>=1.6)', 'MDAnalysis (>=0.16.2)'],
          install_requires=['numpy>=1.6.0', 'MDAnalysis>=0.16.2'],
          provides=['rotcon'],
          license='GPL 2',
          packages=find_packages(exclude=['scripts', 'rotcon/data']),
          package_data={'rotcon': ['data/*.pdb', 'data/*.dcd', 'data/*.dat']},
          scripts=['scripts/convolve-mtss-rotamers.py', 'scripts/convolve-mtss-rotamers_pre.py'],
          long_description_content_type='text/x-rst',
          long_description=LONG_DESCRIPTION,
          project_urls={
              'Documentation': 'https://www.mdanalysis.org/RotamerConvolveMD',
              'Source': 'https://github.com/MDAnalysis/RotamerConvolveMD',
              'Issue Tracker': 'https://github.com/MDAnalysis/RotamerConvolveMD/issues/',
              'Mailing list': 'https://groups.google.com/group/mdnalysis-discussion',
          },
          keywords="science proteins",
          classifiers=['Development Status :: 7 - Inactive',
                       'Environment :: Console',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
                       'Programming Language :: Python',
                       'Programming Language :: Python :: 2',
                       'Programming Language :: Python :: 2.7',
                       'Programming Language :: Python :: 3',
                       'Programming Language :: Python :: 3.4',
                       'Programming Language :: Python :: 3.5',
                       'Programming Language :: Python :: 3.6',
                       'Topic :: Scientific/Engineering',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',
                       'Topic :: Scientific/Engineering :: Chemistry',
    ],
          zip_safe=True,
          )
