
## Run binning tools
<!-- #endregion -->

### metabat2

We first generate a depth file from the bam file using jgi_summarize_bam_contig_depths script from metabat2. This depth file will be used also with maxbin2.    
```bash

jgi_summarize_bam_contig_depths --outputDepth depth_Kickstart.txt alignments_bwa/Kickstart.bam 
```

Now we can run metabat2: 

```bash

metabat2 --inFile Kickstart.megahit/R1.contigs.fa --abdFile depth_Kickstart.txt --outFile metabat2/metabat2 --numThreads 12 --seed  1 

```


### maxbin2

We use the depth file produced by `jgi_summarize_bam_contig_depths`

```bash

mkdir -p maxbin2
run_MaxBin.pl -contig Kickstart.megahit/R1.contigs.fa -abund   depth_Kickstart.txt -thread 12  -out maxbin2/maxbin2

```

### concoct

Then we can also run concoct with the folowing commands:

```bash

mkdir -p concoct/

cut_up_fasta.py Kickstart.megahit/R1.contigs.fa --chunk_size 10000 --overlap_size 0 --merge_last --bedfile concoct/contigs_10K.bed > concoct/contigs_10K.fa

concoct_coverage_table.py concoct/contigs_10K.bed alignments_bwa/Kickstart.bam  > concoct/coverage_table.tsv

concoct --composition_file concoct/contigs_10K.fa --coverage_file concoct/coverage_table.tsv --basename concoct/bins --threads 12

merge_cutup_clustering.py concoct/bins_clustering_gt1000.csv > concoct/clustering_merge.csv

mkdir -p concoct/bins

extract_fasta_bins.py Kickstart.megahit/R1.contigs.fa concoct/clustering_merge.csv --output_path concoct/bins
```

### SemiBin2

We can launch semibin2 as well with its `single_easy_bin` command. 

```{note}
This take some time so it can be skipped.
```

```bash

SemiBin2 single_easy_bin -i Kickstart.megahit/R1.contigs.fa -b alignments_bwa/Kickstart.bam -o semibin2/ -p 12 

```

