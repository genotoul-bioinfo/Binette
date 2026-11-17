
## Run Binette

Binette will use the previously computed bins to refine and improve them, generating a new set of higher-quality bins.

To run Binette, use the following command:

```{include} snippets/05_binette.sh
:code: bash
```

```{admonition} CheckM2 database
:class: tip

This command assumes the CheckM2 database has been downloaded with `checkm2 database --download` as described in the [environment setup](./set_environment.md). If Binette fails with `DIAMOND database not found`, give the database file explicitly with `--checkm2_db /path/to/CheckM2_database/uniref100.KO.1.dmnd`.
```

```{admonition} ⌛ Expected Time
:class: note

This process should take around 5 minutes to complete.
```


Once Binette completes, the `binette_output` directory should have the following structure:

```text
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
├── final_contig_to_bin.tsv
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

```{note}
Some steps of this tutorial (assembly and binning) are not fully deterministic and depend on tool versions, so the number of bins you obtain may differ slightly from the example above.
```

### Key Output Files:

- **`final_bins/`**: Contains the refined bins in FASTA format.
- **`final_bins_quality_reports.tsv`**: A summary report containing CheckM2 metrics for the final bin selection.
- **`final_contig_to_bin.tsv`**: A table assigning each contig to its final bin.
- **`input_bins_quality_reports/`**: Quality reports for each of the input bin sets from MaxBin2, MetaBAT2, CONCOCT, and SemiBin2.

### Next Steps

In the next section, we will use `final_bins_quality_reports.tsv` along with the reports from `binette_output/input_bins_quality_reports` to visualize Binette's bins and compare them with the initial bin sets.

