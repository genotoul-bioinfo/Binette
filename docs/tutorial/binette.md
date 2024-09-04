
## Run Binette

Binette will use the previously computed bins to refine and improve them, generating a new set of higher-quality bins.

To run Binette, use the following command:

```bash
binette --bin_dirs maxbin2/ metabat2/ semibin2/output_bins/ concoct/bins/ \
        -c Kickstart.megahit/R1.contigs.fa \
        --verbose -t 12 -o binette_results
```

```{admonition} ⌛ Expected Time
:class: note
:class: dropdown

This process should talke around 9 minutes to complete.
```


Once Binette completes, the `binette_results` directory should have the following structure:

```
binette_results/
├── final_bins
│   ├── bin_13475.fa
│   ├── bin_17075.fa
│   ├── bin_19689.fa
│   ├── bin_21248.fa
│   ├── bin_31703.fa
│   ├── bin_33569.fa
│   ├── bin_39350.fa
│   ├── bin_39427.fa
│   ├── bin_39558.fa
│   ├── bin_44137.fa
│   ├── bin_46775.fa
│   ├── bin_47060.fa
│   ├── bin_47177.fa
│   ├── bin_47926.fa
│   └── bin_51082.fa
├── final_bins_quality_reports.tsv 
├── input_bins_quality_reports
│   ├── input_bins_1.concoct_bins.tsv
│   ├── input_bins_2.maxbin2.tsv
│   ├── input_bins_3.metabat2.tsv
│   └── input_bins_4.semibin2_output_bins.tsv
└── temporary_files
    ├── assembly_proteins.faa
    ├── diamond_result.log
    └── diamond_result.tsv
```

### Key Output Files:

- **`final_bins/`**: Contains the refined bins in FASTA format.
- **`final_bins_quality_reports.tsv`**: A summary report containing CheckM2 metrics for the final bin selection.
- **`input_bins_quality_reports/`**: Quality reports for each of the input bin sets from MaxBin2, MetaBAT2, CONCOCT, and SemiBin2.

### Next Steps

In the next section, we will use `final_bins_quality_reports.tsv` along with the reports from `binette_results/input_bins_quality_reports` to visualize Binette's bins and compare them with the initial bin sets.

