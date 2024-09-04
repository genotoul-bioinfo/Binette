## Assemble the Reads

We will use **MEGAHIT** to assemble the reads from our dataset. Run the following command:

```{code-block} bash
megahit -1 coal-metagenomics/Kickstart_1.fastq.gz \
        -2 coal-metagenomics/Kickstart_2.fastq.gz \
        --out-dir Kickstart.megahit --out-prefix R1 --num-cpu-threads 12
```

:::{admonition} ⌛ Expected Time
:class: note
:class: dropdown

This process takes approximately 28 minutes to complete.
:::



```{admonition} Assembly tips
:class: tip

Here are some general tips that might help improve your assembly results, depending on your data:

- **Read Cleaning:** If your reads have low-quality bases or adapters, consider cleaning them with a tool like `sickle`. It can boost the overall quality of your assembly.

- **Use SPAdes rather than MEGAHIT:** SPAdes generally performs better than MEGAHIT but takes longer and requires more memory.

- **Quality Check:** Tools like `metaQUAST` are handy for checking your assembly’s quality. It’s a good way to ensure your results are solid before moving on.

- **Assembly Filtering:** After assembling, it’s often a good idea to filter out small or low-coverage contigs. 


These steps aren’t mandatory, and since this tutorial focuses on binning refinement with Binette, we’ll skip them.

```




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
:class: dropdown

This process takes approximately 12 minutes to complete.
:::

```{admonition} Read alignment strategy
:class: tip

If you have multiple samples and assemble them separately, cross-aligning the samples can significantly improve binning. Align each sample to all assemblies and use the resulting BAM files in binning. This approach gives the binning tools more coverage variation, which can be beneficial. However, keep in mind that this process can be resource-intensive, especially with many samples. 

If you did a cross-assembly with your samples, make sure to map the reads separately for each one, generating as many BAM files as you have samples, to help the binning tool. 🚀

```