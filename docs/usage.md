
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

In both formats, the `--contigs` argument should specify a FASTA file containing all the contigs found in the bins. Typically, this file would be the assembly FASTA file used to generate the bins. In these examples the `assembly.fasta` file should contain at least the five contigs mentioned in the `contig2bin_tables` files or in the bin fasta files: `contig_1`, `contig_8`, `contig_15`, `contig_9`, and `contig_10`.



### Providing Precomputed Protein Sequences

You can provide protein sequences in FASTA format to Binette using the `--proteins` argument. The sequence identifiers must follow the Prodigal convention: `<contigID>_<GeneID>`. This naming format ensures proper mapping of each gene to its contig.  

By using this option, the gene prediction step is skipped.  

```{note}
When using precomputed protein sequences, the `coding_density` column in the output reports will be empty, as this metric requires gene coordinates that are only available when genes are freshly predicted.
```

#### Example  
If your contig is named `contig_A`, the gene identifiers should follow this pattern:  
- `contig_A_1`  
- `contig_A_2`  
- `contig_A_3`  


## Outputs

Binette results are stored in the `results` directory. You can specify a different directory using the `--outdir` option.

In this directory you will find:
- `final_bins_quality_reports.tsv`: This is a TSV (tab-separated values) file containing quality information about the final selected bins.
- `final_bins/`: This directory stores all the selected bins in fasta format. Can be skipped with `--no-write-fasta-bins`.
- `final_contig_to_bin.tsv`: A headerless TSV file mapping each contig to its assigned bin. This format is much lighter than the fasta output to describe the final Binette bins.
- `input_bins_quality_reports/`: A directory storing quality reports for the input bin sets, with files following the same structure as `final_bins_quality_reports.tsv`.
- `temporary_files/`: This directory contains intermediate files. If you choose to use the `--resume` option, Binette will utilize files in this directory to prevent the recomputation of time-consuming steps.


The `final_bins_quality_reports.tsv` file contains the following columns:
| Column Name        | Description                                                                                                                                    |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| **name**           | The unique name of the bin.                                                                                                                    |
| **origin**         | Indicates the source of the bin: either an original bin set (e.g., `B`) or `binette` for intermediate bins.                                    |
| **is\_original**   | Boolean flag indicating if the bin is an original bin (`True`) or an intermediate bin (`False`).                                               |
| **original\_name** | The name of the original bin from which this bin was derived.                                                                                  |
| **completeness**   | The completeness of the bin, determined by CheckM2.                                                                                            |
| **contamination**  | The contamination of the bin, determined by CheckM2.                                                                                           |
| **checkm2\_model** | The CheckM2 model used for quality prediction: `Gradient Boost (General Model)` or `Neural Network (Specific Model)`.|
| **score**          | Computed score: `completeness - contamination * weight`. The contamination weight can be customized using the `--contamination_weight` option. |
| **size**           | Total size of the bin in nucleotides.                                                                                                          |
| **N50**            | The N50 of the bin, representing the length for which 50% of the total nucleotides are in contigs of that length or longer.                    |
| **coding\_density** | The percentage of the bin that codes for proteins (genes length / total bin length × 100). Only computed when genes are freshly identified. Empty when using `--proteins` or `--resume` options. |
| **contig\_count**  | Number of contigs contained within the bin.                                                                                                    |
   


