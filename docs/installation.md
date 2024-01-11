
# Installation

## With Bioconda

Binette can be esailly installed with conda 

```bash
conda create -c bioconda -c defaults -c conda-forge -n binette binette
conda activate binette
```

Binette should be able to run :

```
binette -h
```


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

Finally install Binette with pip

```
pip install .
```

Binette should be able to run :

```
binette -h
```


## Downloading the CheckM2 database

Before using Binette, it is necessary to download the CheckM2 database:

```bash
checkm2 database --download --path <checkm2/database/>
```

Make sure to replace `<checkm2/database/>` with the desired path where you want to store the CheckM2 database.
