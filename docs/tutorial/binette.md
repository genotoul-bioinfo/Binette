
## Run Binette

Binette will use the previously computed bins to refine and improve them, generating a new set of higher-quality bins.

To run Binette, use the following command:

```{include} snippets/05_binette.sh
:code: bash
```

```{admonition} ⌛ Expected Time
:class: note

This process should talke around 5 minutes to complete.
```


Once Binette completes, the `binette_results` directory should have the following structure:

```
binette_output/
├── final_bins
│   ├── binette_bin10.fa
│   ├── binette_bin11.fa
│   ├── binette_bin12.fa
│   ├── binette_bin13.fa
│   ├── binette_bin14.fa
│   ├── binette_bin15.fa
│   ├── binette_bin1.fa
│   ├── binette_bin2.fa
│   ├── binette_bin3.fa
│   ├── binette_bin4.fa
│   ├── binette_bin5.fa
│   ├── binette_bin6.fa
│   ├── binette_bin7.fa
│   ├── binette_bin8.fa
│   └── binette_bin9.fa
├── final_bins_quality_reports.tsv
├── input_bins_quality_reports
│   ├── input_bins_1.concoct_bins.tsv
│   ├── input_bins_2.maxbin2_bins.tsv
│   ├── input_bins_3.metabat2_bins.tsv
│   └── input_bins_4.semibin2_output_output_bins.tsv
└── temporary_files
    ├── assembly_proteins.faa.gz
    ├── diamond_result.log
    └── diamond_result.tsv.gz

```

### Key Output Files:

- **`final_bins/`**: Contains the refined bins in FASTA format.
- **`final_bins_quality_reports.tsv`**: A summary report containing CheckM2 metrics for the final bin selection.
- **`input_bins_quality_reports/`**: Quality reports for each of the input bin sets from MaxBin2, MetaBAT2, CONCOCT, and SemiBin2.

### Next Steps

In the next section, we will use `final_bins_quality_reports.tsv` along with the reports from `binette_results/input_bins_quality_reports` to visualize Binette's bins and compare them with the initial bin sets.

