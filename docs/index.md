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


Binette is a fast and accurate binning refinement tool to constructs high quality MAGs from the output of multiple binning tools.

From the input bin sets, Binette constructs new hybrid bins. A bin can be seen as a set of contigs. When at least two bins overlap, meaning they share at least one contig, Binette utilizes basic set operations to create new bins.
- Intersection bin: This bin consists of the contigs that are shared by the overlapping bins. 
- Difference bin: This bin contains the contigs that are exclusively found in one bin and not present in the others.
- Union bin: The union bin includes all the contigs contained within the overlapping bins

It then uses CheckM2 to assess bins quality to finally select the best bins possible.

Binette is inspired from the metaWRAP bin-refinement tool but it effectively solves all the problems from that very tool. 
- Enhanced Speed: Binette significantly improves the speed of the refinement process. It achieves this by launching the initial steps of CheckM2, such as Prodigal and Diamond runs, only once on all contigs. These intermediate results are then utilized to assess the quality of any given bin, eliminating redundant computations and accelerating the refinement process.
- No Limit on Input Bin Sets: Unlike its predecessor, Binette is not constrained by the number of input bin sets. It can handle and process multiple bin sets simultaneously.




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

