
## Align the Reads to the Assembly

Binning tools rely on coverage information, among other criteria, to evaluate each contig. 

To obtain this coverage data, we first need to map the reads back to the assembly.

```{code-block} bash
# Create a directory for the alignments
mkdir -p alignments_bwa/

# Index the contigs file using BWA-MEM2
bwa-mem2 index Kickstart.megahit/R1.contigs.fa -p Kickstart.megahit/R1.contigs.fa

# Map reads back to the assembly, convert to BAM format, and sort
bwa-mem2 mem -t 12 Kickstart.megahit/R1.contigs.fa coal-metagenomics/Kickstart_*.fastq.gz | \
samtools view -@ 12 -bS - | \
samtools sort -@ 12 - -o alignments_bwa/Kickstart.bam

# Index the BAM file
samtools index alignments_bwa/Kickstart.bam
```


:::{admonition} ⌛ Expected Time
:class: note

This process takes approximately 12 minutes to complete.
:::

```{admonition}
:class: tip

If you have multiple samples and assemble them separately, cross-aligning the samples can significantly improve binning. Align each sample to all assemblies and use the resulting BAM files in binning. This approach gives the binning tools more coverage variation, which can be beneficial. However, keep in mind that this process can be resource-intensive, especially with many samples. 

If you did a cross-assembly with your samples, make sure to map the reads separately for each one, generating as many BAM files as you have samples, to help the binning tool. 🚀

```


## Run Binning Tools

Let's use different binning tools to group the contigs into bins, which we'll refine in the next section with Binette.

### MetaBAT2

First, generate a depth file from the BAM file using the `jgi_summarize_bam_contig_depths` script from MetaBAT2. This depth file will also be used by MaxBin2. 

```bash
jgi_summarize_bam_contig_depths --outputDepth depth_Kickstart.txt alignments_bwa/Kickstart.bam
```

Now, run MetaBAT2 with the generated depth file:

```bash
metabat2 --inFile Kickstart.megahit/R1.contigs.fa --abdFile depth_Kickstart.txt --outFile metabat2/metabat2 --numThreads 12 --seed 1
```

### MaxBin2

We will use the same depth file produced by `jgi_summarize_bam_contig_depths` for MetaBAT2:

```bash
mkdir -p maxbin2
run_MaxBin.pl -contig Kickstart.megahit/R1.contigs.fa \
                -abund depth_Kickstart.txt -thread 12 -out maxbin2/maxbin2
```

### CONCOCT

To run CONCOCT, follow these steps:

1. **Cut up the FASTA file** into chunks for processing:

```bash
mkdir -p concoct/

cut_up_fasta.py Kickstart.megahit/R1.contigs.fa --chunk_size 10000 \
                --overlap_size 0 --merge_last \
                --bedfile concoct/contigs_10K.bed > concoct/contigs_10K.fa
```

2. **Generate the coverage table** from the BAM file:

```bash
concoct_coverage_table.py concoct/contigs_10K.bed alignments_bwa/Kickstart.bam > concoct/coverage_table.tsv
```

3. **Run CONCOCT** with the composition and coverage files:

```bash
concoct --composition_file concoct/contigs_10K.fa \
        --coverage_file concoct/coverage_table.tsv \
        --basename concoct/bins --threads 12
```

4. **Merge the clustering results** and extract bins:

```bash
merge_cutup_clustering.py concoct/bins_clustering_gt1000.csv > concoct/clustering_merge.csv

mkdir -p concoct/bins

extract_fasta_bins.py Kickstart.megahit/R1.contigs.fa concoct/clustering_merge.csv --output_path concoct/bins
```

### SemiBin2

You can also run SemiBin2 with its `single_easy_bin` command:

```{admonition} ⏳ Time Note
:class: note

This process can take some time, so it may be skipped.
```

```bash
SemiBin2 single_easy_bin -i Kickstart.megahit/R1.contigs.fa \
                            -b alignments_bwa/Kickstart.bam \
                            -o semibin2/ -p 12
```
