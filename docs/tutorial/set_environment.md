## Set Up the Tutorial Environment

To get started, we'll download the necessary tools and set them up in a dedicated Conda environment.

### Create a Conda Environment

First, let's create a new Conda environment specifically for this tutorial:

```{code-block} bash
mamba env create -f binette_tutorial_env.yaml -n binette_tutorial
```

This command will create a Conda environment named `binette_tuto` using the environment file `binette_tutorial_env.yaml`.

Below is the content of the `binette_tutorial_env.yaml` file:

```{include} binette_tutorial_env.yaml
:code: yaml
```


### Activate the Environment

After the environment is created, activate it by running:

```{code-block} bash
conda activate binette_tuto
```

