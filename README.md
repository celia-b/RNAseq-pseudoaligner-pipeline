# RNAseq-pseudoaligner-pipeline

This pipeline allows the user to process raw RNA-seq data in an automatized fashion, producing a quantification that is ready for differential expression analysis. 

**Basic structure of the pipeline**

Here is a breakdown of the pipeline structure and steps:

1. Quality check of the raw reads & subsequent trimming. This program is not able to read the content of the quality report, so it performs trimming regardless, as trimming good quality reads leads to no significant change in them. _Note: the trimming step doesn't currently work. Must be run separately_
2. Quality check of the processed reads.
3. Transcriptome download and indexing by kallisto.
4. Creation of a transcript-to-gene translation table. We want to be able to connect transcript names with their genes of origin, and this is what this table does.
5. Quantification by kallisto.
6. Augmentation of quantification table with gene names (as provided by the translation table in step 4).
7. Merging of all samples into a single table, for further processing in DESeq2.

Steps 1-4 can be manually cancelled if they have been run independently. For that, open the pipeline script and read the section with header _"README: User defined parameters"_.
Furthermore, a DESeq2-based R notebook is provided with the main steps and code needed for basic DE analysis of the produced results. For more information on the package, see <http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#se>

**Underlying tools**

It is built around the following tools:

- fastqc & multiqc: for quality control
- trimommatic: for read trimming
- kallisto: for rapid pseudoalignment of reads and quantification

**Parallelization**

The pipeline runs each sample in parallel, so runtime is relatively fast when compared to a sequential run. Because of memory limitations, no more than 20 samples can be processed simultaneously. Given that typical kallisto run on a 30M read sample takes <10 minutes, any experiment with <=20 samples will also take <10 minutes for kallisto. Extra samples will progressively increase the runtime. Further steps in the pipeline also increase runtime, with 40 samples taking around 25 minutes in total.


**Directory tree**

The user provides a project directory, which must contain a folder called raw_reads/ containing the raw reads of each sample, each within their own folder (see diagram). The output of the pipeline is organized as shown in the diagram below, with folders for the QC, trimming and quantification with kallisto. Some files like the transcriptome, the index, the translation table, the final merged samples table and the logs live directly in the project folder.

_Note: it has been made possible for the pipeline to take trimmed reads from the raw_reads/ folder, only because this was the original disposition of the files. However, when the pipeline runs the trimming itself, files go into the trimmed_reads/ folder. For future projects and if trimming is performed outside of the pipeline, it is suggested to put those files in the trimmed_reads/ folder, for obvious reasons._

RNAseq-project/  
├── raw_reads/  
│   ├── S100/  
│   │   ├── S100_1.fq.gz
│   │   ├── S100_2.fq.gz
│   │   ├── S100_forward_paired.fq.gz (optional, can also be in trimmed_reads/)
│   │   ├── S100_forward_unpaired.fq.gz (optional, can also be in trimmed_reads/)
│   │   ├── S100_reverse_paired.fq.gz (optional, can also be in trimmed_reads/)
│   │   └── S100_reverse_unpaired.fq.gz (optional, can also be in trimmed_reads/)
│   └── ...
├── trimmed_reads/   
│   ├── S100/
│   │   ├── S100_forward_paired.fq.gz (optional, can also be in raw_reads/)
│   │   ├── S100_forward_unpaired.fq.gz (optional, can also be in raw_reads/)
│   │   ├── S100_reverse_paired.fq.gz (optional, can also be in raw_reads/)
│   │   └── S100_reverse_unpaired.fq.gz (optional, can also be in raw_reads/)
│   └── ...
├── QC/
│   ├── raw_reads/
│   │   ├── S100_1_fastqc.html
│   │   ├── S100_2_fastqc.zip
│   │   └── ...
│   └── trimmed_reads/
│       ├── S100_forward_paired_fastqc.html
│       ├── S100_forward_paired_fastqc.zip
│       └── ...
├── transcriptome_file.fa.gz
├── indexed_transcriptome
├── translation_table.json
├── quantification_kallisto/
│   ├── S100/
│   │   ├── abundance_aug.tsv
│   │   ├── abundance.tsv
│   │   ├── abundance.h5
│   │   └── run_info.json
│   └── ...
├── merged_samples.tsv
├── kallisto.err
└── kallisto.log


