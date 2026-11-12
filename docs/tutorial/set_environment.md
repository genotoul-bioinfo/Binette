## Set Up the Tutorial Environment

To get started, we'll download the necessary tools and set them up in a dedicated Conda environment.

### Create a Conda Environment

First, download the environment file `binette_tutorial_env.yaml`:

```{code-block} bash
wget https://raw.githubusercontent.com/genotoul-bioinfo/Binette/main/docs/tutorial/binette_tutorial_env.yaml
```

Then create a new Conda environment specifically for this tutorial:

```{code-block} bash
conda env create -f binette_tutorial_env.yaml
```

This command will create a Conda environment named `binette_tutorial` using the environment file `binette_tutorial_env.yaml`.

Below is the content of the `binette_tutorial_env.yaml` file:

```{include} binette_tutorial_env.yaml
:code: yaml
```

### Activate the Environment

After the environment is created, activate it by running:

```{code-block} bash
conda activate binette_tutorial
```

### Download the CheckM2 Database

Binette uses the CheckM2 database to assess bin quality. With the environment activated, download it (~3 GB) with:

```{code-block} bash
checkm2 database --download --path checkm2_database/
```

This downloads the database and registers its location in the environment, so Binette will find it automatically.

:::{admonition} Already have the CheckM2 database?
:class: tip

If you have already downloaded the CheckM2 database elsewhere, you don't need to download it again. Either register it in the environment:

```{code-block} bash
checkm2 database --setdblocation /path/to/CheckM2_database/uniref100.KO.1.dmnd
```

or give it directly to Binette with the `--checkm2_db` option when running it in the [Run Binette](./binette.md) section.
:::
