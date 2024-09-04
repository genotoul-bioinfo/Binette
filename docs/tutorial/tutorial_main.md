
# Tutorial 

In this tutorial, we'll walk through a practical example of how to use Binette with real data. We'll start by downloading metagenomics reads and then assemble these reads into contigs. Next, we'll use different binning tools to group the contigs into bins. Finally, we'll use Binette to refine these bins.

```{mermaid}

--- 
title: "Tutorial Overview:"
align: center

config:
  look: handDrawn
  theme: neutral
---


graph TD

    i[metagenomics reads] --> B[assembly]


    B --> metabat2 --> r[binette]
    B --> maxbin2 --> r
    B --> concoct --> r
    B --> semibin2 --> r
    r --> f[final bins]
    
        subgraph Binning
            metabat2:::binning
            maxbin2:::binning
            concoct:::binning
            semibin2:::binning
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
analyse_binette_result.myst
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
