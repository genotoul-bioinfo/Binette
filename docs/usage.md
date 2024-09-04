
# Usage 

## Input Formats

Binette supports two input formats for bin sets: 

1. **Contig2bin Tables:** You can provide bin sets using contig2bin tables, which establish the relationship between each contig and its corresponding bin. In this format, you need to specify the `--contig2bin_tables` argument. 

For example, consider the following two `contig2bin_tables`:

- `bin_set1.tsv`:

    ```
    contig_1   binA
    contig_8   binA
    contig_15  binB
    contig_9   binC
    ```
    
- `bin_set2.tsv`:

    ```
    contig_1   bin.0
    contig_8   bin.0
    contig_15  bin.1
    contig_9   bin.2
    contig_10  bin.0
    ```
    
    The `binette` command to process this input would be:
    
    ```bash
    binette --contig2bin_tables bin_set1.tsv bin_set2.tsv --contigs assembly.fasta
    ```

2. **Bin Directories:** Alternatively, you can use bin directories, where each bin is represented by a separate FASTA file. For this format, you need to provide the `--bin_dirs` argument. Here's an example of two bin directories:

    ```
    bin_set1/
    ├── binA.fa: contains sequences of contig_1, contig_8
    ├── binB.fa: contains sequences of contig_15
    └── binC.fa: contains sequences of contig_9
    ```
    
    ```
    bin_set2/
    ├── binA.fa: contains sequences of contig_1, contig_8, contig_10
    ├── binB.fa: contains sequences of contig_15
    └── binC.fa: contains sequences of contig_9
    ```
    
    The `binette` command to process this input would be:
    
    ```bash
    binette --bin_dirs bin_set1 bin_set2 --contigs assembly.fasta
    ```

In both formats, the `--contigs` argument should specify a FASTA file containing all the contigs found in the bins. Typically, this file would be the assembly FASTA file used to generate the bins. In these exemple the `assembly.fasta` file should contain at least the five contigs mentioned in the `contig2bin_tables` files or in the bin fasta files: `contig_1`, `contig_8`, `contig_15`, `contig_9`, and `contig_10`.

## Outputs

Binette results are stored in the `results` directory. You can specify a different directory using the `--outdir` option.

In this directory you will find:
- `final_bins_quality_reports.tsv`: This is a TSV (tab-separated values) file containing quality information about the final selected bins.
- `final_bins/`: This directory stores all the selected bins in fasta format.
- `input_bins_quality_reports/`: A directory storing quality reports for the input bin sets, with files following the same structure as `final_bins_quality_reports.tsv`.
- `temporary_files/`: This directory contains intermediate files. If you choose to use the `--resume` option, Binette will utilize files in this directory to prevent the recomputation of time-consuming steps.


The `final_bins_quality_reports.tsv` file contains the following columns:
| Column Name         | Description                                                                                                  |
|---------------------|--------------------------------------------------------------------------------------------------------------|
| **bin_id**          | This column displays the unique ID of the bin.                                                             |
| **origin**          | Indicates the source or origin of the bin, specifying from which bin set it originates or the intermediate set operation that created it. |
| **name**            | The name of the bin.                                                                                        |
| **completeness**    | The completeness of the bin, determined by CheckM2.                                                         |
| **contamination**   | The contamination of the bin, determined by CheckM2.                                                       |
| **score**           | This column displays the computed score, which is calculated as: `completeness - contamination * weight`. You can customize the contamination weight using the `--contamination_weight` option. |
| **size**            | Represents the size of the bin in nucleotides.                                                              |
| **N50**             | Displays the N50 of the bin.                                                                                |
| **contig_count**    | The number of contigs contained within the bin.          