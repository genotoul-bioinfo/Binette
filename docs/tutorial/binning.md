## Run Binning Tools

Let's use different binning tools to group the contigs into bins, which we'll refine in the next section with Binette.

### MetaBAT2

```{include} snippets/04a_binning_metabat2.sh
:code: bash
```

### MaxBin2

```{include} snippets/04b_binning_maxbin2.sh
:code: bash
```

### CONCOCT

```{include} snippets/04c_binning_concoct.sh
:code: bash
```

```{note}
CONCOCT may print a `SyntaxWarning: invalid escape sequence` and many `No consensus cluster for contig ...` messages. These are expected and can be safely ignored.
```

### SemiBin2

```{include} snippets/04d_binning_semibin2.sh
:code: bash
```

```{admonition} ⌛ Expected Time
:class: note

With 12 threads, MetaBAT2 and MaxBin2 take a couple of minutes, CONCOCT around 15 minutes and SemiBin2 around 1 hour.
```
