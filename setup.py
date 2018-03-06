from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

if sys.version_info < (3,4) and sys.version_info < (2.7):
   sys.exit("Error: You are using Python "+str(sys.version_info)+"; Please use Python  3.4 or 2.7 or better\n")

this_folder = path.abspath(path.dirname(__file__))
with open(path.join(this_folder,'README.md'),encoding='utf-8') as inf:
  long_description = inf.read()

setup(
  name='GSVA',
  version='1.0.5',
  description='Python CLI and module for running the GSVA R bioconductor package with Python Pandas inputs and outputs.',
  long_description=long_description,
  url='https://github.com/jason-weirather/GSVA',
  author='Jason L Weirather',
  author_email='jason.weirather@gmail.com',
  license='Apache License, Version 2.0',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Apache Software License'
  ],
  keywords='bioinformatics, R, enrichment, GSVA, ssGSEA, GSEA, bioconductor',
  packages=['GSVA'
           ],
  package_data={'GSVA':['*.r']},
  install_requires=['pandas'],
  entry_points = {
    'console_scripts':['GSVA=GSVA:__cli']
  }
)

