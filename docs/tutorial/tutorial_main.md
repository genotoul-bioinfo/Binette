
# Tutorial 

In this tutorial, we'll walk through a practical example of how to use Binette with real data.

  1. We'll start by downloading metagenomics reads and then assemble these reads into contigs.
  2. Next, we'll use different binning tools to group the contigs into bins.
  3. Finally, we'll use Binette to refine these bins.

```{admonition} Computational resources
:class: note

The commands in this tutorial use 12 threads (`-t 12`, `--num-cpu-threads 12`, ...): adjust this value to your machine.
The full analysis requires around 6 GB of disk space, plus ~3 GB for the CheckM2 database, and takes about 2.5 hours in total with 12 threads.
```

```{mermaid}

--- 
title: "Tutorial Overview:"
align: right

config:
  look: handDrawn
  theme: neutral
---


graph TD

    i[Metagenomics reads] --> B[Assembly]


    B --> MetaBAT2 --> r[Binette]
    B --> MaxBin2 --> r
    B --> CONCOCT --> r
    B --> SemiBin2 --> r
    r --> f[final bins]
    
        subgraph Binning
            MetaBAT2:::binning
            MaxBin2:::binning
            CONCOCT:::binning
            SemiBin2:::binning
        end

            
        classDef binning fill:#d4ae40


```


```{toctree}
:caption: 'Tutorial steps'
:maxdepth: 1

set_environment
get_dataset
assembly
binning
binette
analyse_binette_result.ipynb
```


<!-- 
```{include} ./set_env_and_get_data.md
```



```{include} ./assembly.md
```

```{include} ./binning.md
```

```{include} ./binette.md
```

```{include} ./analyse_binette_result.ipynb
``` -->


<!-- ### Compare binette results with intial bins



## Compare Binette and Das Tool 


### Run Das Tool 


### Compare binette results with intial bins


 -->
