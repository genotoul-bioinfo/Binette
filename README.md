# Overview 

Binette is a fast and accurate binning refinement tool that constructs high quality MAGs from the output of multiple binning tools.

From the input bin sets, Binette constructs new hybrid bins. A bin can be seen as a set of contigs. Based on this property when a least two bins overlap (they share at least one contig), Binette creates new bins using basic set operations:
- Intersection bin: contigs shared by the overlapping bins.
- Difference bin: contigs that are found only in one bin and not in the other ones.
- Union bin : all contigs contained in the overlapping bins.

It then uses checkm2 to assess bins quality to finally select the best bins possible.

Binette is inspired from the metaWRAP bin-refinement tool but it effectively solves all the problems from that very tool. 
- It is much faster as it launches the first steps of checkm2 (prodigal and diamond runs) only once on all contigs and then uses this intermediate results to assess quality for any bin.
- It is not limited in the number of input bin sets.
- It selects the best bins in a more accurate and elegant manner.
- It is easier to use.

# Installing

You can install binette directly from the source code or build and run it from within a Singularity container.

## Installing from a conda environnement

Clone this repository: 
```
git clone https://github.com/genotoul-bioinfo/Binette
cd Binette
```

Then create a Conda environment using the `binette_dev.yaml` file:
```
conda env create -n binette -f binette_dev.yaml
conda activate binette 
```

Binette need checkm2 to be fully installed with pip.

Follow Chekm2 installation instruction:

```
git clone --recursive https://github.com/chklovski/checkm2.git

pip install checkm2/

```
and download its database

```

checkm2 database --download --path checkm2/database/
```

Finally install binette with pip

```
pip install .
```

Binette should be able to run :

```
binette -h
```


## Running binette using a singularity image

Singularity version 3 or above must be installed. See [here](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps) how to install Singularity >=v3.

Git clone binette repository and build the singularity image. 

```
git clone https://forgemia.inra.fr/jean.mainguy/binette
cd binette
sudo singularity build binette.sif singularity_recipe
```

Then if the build succesfully finished, you should be able to run Binette:

```
singularity exec binette.sif python binette/binette.py -h
```

# Usage 



# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker:

[binette issue tracker](https://github.com/genotoul-bioinfo/Binette/issues)

# Licence

This program is released as an open source software under the terms of [MIT License](https://forgemia.inra.fr/jean.mainguy/binette/-/raw/main/LICENSE).
