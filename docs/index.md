% Binette documentation master file, created by
% sphinx-quickstart on Thu Jan 11 21:13:20 2024.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.


[![Bioconda](https://img.shields.io/conda/vn/bioconda/binette?label=bioconda&logo=anaconda&logoColor=white)](https://bioconda.github.io/recipes/binette/README.html)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/binette?label=bioconda%20downloads)](https://anaconda.org/bioconda/binette)
[![PyPI](https://img.shields.io/pypi/v/binette?logo=pypi&logoColor=white)](https://pypi.org/project/Binette/)
[![Python](https://img.shields.io/pypi/pyversions/binette?logo=python&logoColor=white)](https://pypi.org/project/Binette/)
[![License](https://img.shields.io/github/license/genotoul-bioinfo/Binette)](https://github.com/genotoul-bioinfo/Binette/blob/main/LICENSE)

[![CI](https://img.shields.io/github/actions/workflow/status/genotoul-bioinfo/Binette/binette_ci.yml?branch=main&label=CI&logo=github)](https://github.com/genotoul-bioinfo/Binette/actions/workflows/binette_ci.yml)
[![Coverage](https://genotoul-bioinfo.github.io/Binette/coverage-badge.svg)](https://genotoul-bioinfo.github.io/Binette/)
[![Documentation](https://img.shields.io/readthedocs/binette/latest?logo=readthedocs&logoColor=white)](https://binette.readthedocs.io/en/latest/)

[![JOSS](https://joss.theoj.org/papers/10.21105/joss.06782/status.svg)](https://doi.org/10.21105/joss.06782)
[![Source](https://img.shields.io/badge/source-GitHub-303030?logo=github)](https://github.com/genotoul-bioinfo/Binette)

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

