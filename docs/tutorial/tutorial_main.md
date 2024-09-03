
# Tutorial 

In this tutorial, we'll walk through a practical example of how to use Binette with real data. We'll start by downloading metagenomics reads and then assemble these reads into contigs. Next, we'll use different binning tools to group the contigs. Finally, we'll use Binette to refine these bins and improve our results.

```{mermaid}
---
title: "Tutorial Overview:"
align: center
---

%%{init: {'theme':'default'}}%%

graph LR

  A[Download Metagenomics Reads] --> B
  B[Assemble Reads into Contigs] --> c
          subgraph Pangenome creation
            a:::workflow
            c:::workflow
            g:::workflow
            p:::workflow
            a("annotate") --> c
            c(cluster) --> g(graph)
            g(graph) --> p(partition)
        end


  C[Bin Contigs with Binning Tools] --> D[Refine Bins with Binette]


        
    classDef panrgp fill:#4066d4
    classDef panmodule fill:#d44066
    classDef workflow fill:#d4ae40


```

```{mermaid}

---
title: "Tutorial Overview:"
align: center
---


graph TD

    i[Get Metagenomics Reads] --> B[Assembly & Reads alignment]


    B --> metabat2 --> r[Binette]
    B --> maxbin2 --> r
    B --> concoct --> r
    B --> semibin2 --> r
    
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
