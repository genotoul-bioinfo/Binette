# Overview 

Binette is a fast and accurate binning refinement tool to constructs high quality MAGs from the output of multiple binning tools.

From the input bin sets, Binette constructs new hybrid bins. A bin can be seen as a set of contigs. When at least two bins overlap, meaning they share at least one contig, Binette utilizes basic set operations to create new bins.
- Intersection bin: This bin consists of the contigs that are shared by the overlapping bins. 
- Difference bin: This bin contains the contigs that are exclusively found in one bin and not present in the others.
- Union bin: The union bin includes all the contigs contained within the overlapping bins

It then uses checkm2 to assess bins quality to finally select the best bins possible.

Binette is inspired from the metaWRAP bin-refinement tool but it effectively solves all the problems from that very tool. 
- Enhanced Speed: Binette significantly improves the speed of the refinement process. It achieves this by launching the initial steps of checkm2, such as prodigal and diamond runs, only once on all contigs. These intermediate results are then utilized to assess the quality of any given bin, eliminating redundant computations and accelerating the refinement process.
- No Limit on Input Bin Sets: Unlike its predecessor, Binette is not constrained by the number of input bin sets. It can handle and process multiple bin sets simultaneously, accommodating a broader range of data and experimental setups.
<!-- - Bin selection have been improved. It selects the best bins in a more accurate and elegant manner.
- It is easier to use. -->

# Installation

## From a conda environnement

Clone this repository: 
```
git clone https://github.com/genotoul-bioinfo/Binette
cd Binette
```

Then create a Conda environment using the `binette.yaml` file:
```
conda env create -n binette -f binette.yaml
conda activate binette 
```
<!-- 
Binette need checkm2 to be fully installed with pip.

Follow Chekm2 installation instruction:

You can install it with git: 

```
git clone --recursive https://github.com/chklovski/checkm2.git

pip install checkm2/

```
Or download the archive from github:

```bash
# get the archive
wget https://github.com/chklovski/CheckM2/archive/refs/tags/1.0.2.tar.gz

# decompress
tar -xf 1.0.2.tar.gz
rm 1.0.2.tar.gz

# install
pip install CheckM2-1.0.2/

``` -->

Download checkm2 database

```

checkm2 database --download --path <checkm2/database/>
```

Finally install binette with pip

```
pip install .
```

Binette should be able to run :

```
binette -h
```


<!-- ## Running binette using a singularity image

Singularity version 3 or above must be installed. See [here](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps) how to install Singularity >=v3.

Git clone binette repository and build the singularity image. 

```
git clone https://github.com/genotoul-bioinfo/Binette
cd Binette
sudo singularity build binette.sif singularity_recipe
```

Then if the build succesfully finished, you should be able to run Binette:

```
singularity exec binette.sif binette -h
``` -->

# Usage 



# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker:

[binette issue tracker](https://github.com/genotoul-bioinfo/Binette/issues)

# Licence

This program is released as an open source software under the terms of [MIT License](https://forgemia.inra.fr/jean.mainguy/binette/-/raw/main/LICENSE).
