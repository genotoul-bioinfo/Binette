% Binette documentation master file, created by
% sphinx-quickstart on Thu Jan 11 21:13:20 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.


[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/binette/README.html)  [![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/downloads.svg)](https://anaconda.org/bioconda/binette)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/license.svg)](https://anaconda.org/bioconda/binette) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/binette/badges/version.svg)](https://anaconda.org/bioconda/binette)
[![PyPI version](https://badge.fury.io/py/Binette.svg)](https://badge.fury.io/py/Binette)


[![Test Coverage](https://genotoul-bioinfo.github.io/Binette/coverage-badge.svg)](https://genotoul-bioinfo.github.io/Binette/) 
[![CI Status](https://github.com/genotoul-bioinfo/Binette/actions/workflows/binette_ci.yml/badge.svg)](https://github.com/genotoul-bioinfo/Binette/actions/workflows)
[![Documentation Status](https://readthedocs.org/projects/binette/badge/?version=latest)](https://binette.readthedocs.io/en/latest/?badge=latest)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/genotoul-bioinfo/Binette)

# Binette

**Binette** is a fast and accurate binning refinement tool designed to construct high-quality MAGs from the output of multiple binning tools.  

### How It Works  

From the input bin sets, Binette constructs new **hybrid bins**. A **bin** is a set of contigs. When at least two bins overlap (i.e., share at least one contig), Binette applies fundamental set operations to generate new bins:  

- **Intersection Bin**: Contains contigs that are shared by overlapping bins.  
- **Difference Bin**: Includes contigs that are unique to one bin and not found in others.  
- **Union Bin**: Encompasses all contigs from the overlapping bins.  

Binette then evaluates bin quality using **CheckM2**, ensuring the selection of the best possible bins.  

### Why Binette?  

Binette is inspired by the **metaWRAP bin-refinement tool** but effectively addresses its limitations. Key improvements include:  

- **Enhanced Speed**  
  Binette significantly accelerates the refinement process by running **CheckM2's** initial steps (e.g., Prodigal and Diamond) only once for all contigs. These intermediate results are then reused for bin quality assessment, eliminating redundant computations.  

- **No Limit on Input Bin Sets**  
  Unlike metaWRAP, Binette supports **any number of input bin sets**, allowing seamless processing of multiple binning outputs.  




```{toctree}
:caption: 'Documentation'
:maxdepth: 2

installation
usage
tutorial/tutorial_main
contributing
tests
api/api_ref
```

