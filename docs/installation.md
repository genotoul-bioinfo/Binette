
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


```{tip}
For quicker installation and potential resolution of conflicting dependencies, consider using [Mamba](https://github.com/mamba-org/mamba), an efficient alternative to conda.

```


## Installing from Source Code within a conda environnement

A straightforward method to install Binette from the source code is by utilizing a conda environment that includes all the necessary dependencies.

**1. Clone the Binette Repository**

```bash
git clone https://github.com/genotoul-bioinfo/Binette
cd Binette
```

**2. Installing Dependencies with a Conda Environment File**

Install Binette dependencies listed in the [binette.yaml](https://github.com/genotoul-bioinfo/Binette/blob/main/binette.yaml) file located at the root of the repository, using conda:

```bash
conda env create -n binette -f binette.yaml
conda activate binette
```

**3. Installing Binette**

Finally, install Binette using **pip**:

```bash
pip install .
```

Binette should be able to run :

```bash
binette -h
```


## Downloading the CheckM2 database

Before using Binette, it is necessary to download the CheckM2 database:

```bash
checkm2 database --download --path <checkm2/database/>
```

Make sure to replace `<checkm2/database/>` with the desired path where you want to store the CheckM2 database.
