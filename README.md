# RNAseq-pseudoaligner-pipeline

This pipeline allows the user to process raw RNA-seq reads in an automatized fashion, producing a quantification that is ready for differential expression analysis. 

<u>Basic structure of the pipeline</u>
Here is a breakdown of the pipeline structure and steps:

1. Quality check of the raw reads & subsequent trimming. This program is not able to read the content of the quality report, so it performs trimming regardless, as trimming good quality reads leads to no significant change in them. 
2. Quality check of the processed reads.
3. Transcriptome download and indexing by kallisto.
4. Creation of a translation table. We want to be able to connect transcript names with their genes of origin, this is what this table does.
5. Quantification by kallisto.
6. Augmentation of quantification table with gene names (as provided by the translation table in 4).
7. Merging of samples in a single table, for further processing in DESeq2.

Steps 1-4 can be manually cancelled if they have been run independently. For that, open the pipeline script and read the section with header _"README: User defined parameters"_.
Furthermore, a DESeq2-based R notebook is provided with the main steps and code needed for basic DE analysis of the produced results.

<u>Underlying tools>
It is built around the following tools:

- fastqc & multiqc: for quality control
- trimommatic: for read trimming
- kallisto: for rapid pseudoalignment of reads and quantification

<u>Parallelization</u>
It does so parallelizing the samples, so runtime is relatively fast when compared to a sequential run. Because of memory limitations, no more than 20 samples can be processed simultaneously. Given that typical kallisto run on a 30M read sample takes <10 minutes, any experiment with <=20 samples will also take <10 minutes for kallisto. Extra samples will progressively increase the runtime. Further steps in the pipeline also increase runtime, with 40 samples taking around 25 minutes in total.






