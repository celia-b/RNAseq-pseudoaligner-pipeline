#!/bin/bash
### ----------------------------------- PBS SETTINGS ----------------------------------- ###
### Account information
#PBS -W group_list=dtu_00001 -A dtu_00001
### Job name
#PBS -N kallisto
### Output files
#PBS -e kallisto.err
#PBS -o kallisto.log
### Number of nodes
#PBS -l nodes=1:ppn=20
### Memory
#PBS -l mem=160gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> 
#PBS -l walltime=06:00:00
### ------------------------------------------------------------------------------------ ###




### -------------------------- README + USER-DEFINED PARAMETERS ------------------------- ###
## PROJECT FOLDER
# Please provide a project folder. This folder must contain as many subfolders 
# as there are samples, called S+number (e.g. S1, S2, S100...).
# Inside of these folders, the raw fastq files must be located. 
#
# Please state the path to the project folder:
project_folder="/home/projects/dtu_00001/data/TORSK-CELIA/RNAseq-project/"

## PIPELINE FOLDER
# Please specify the folder where this pipeline & all supplementary scripts are located
pipeline_folder="/home/projects/dtu_00001/data/TORSK-CELIA/kallisto_pipeline/"
#
# With that we can find all necessary scripts - no need to modify this
make_translation=$pipeline_folder"make_translation.py"
add_gene_name=$pipeline_folder"add_gene_name.py"
merge_samples=$pipeline_folder"merge_samples.py"

## STEPS
# This pipeline goes through a number of steps, from quality control to quantification.
# Some of these steps might have previously been completed, so the user can ask for them
# to be skipped in order to reduce runtime. Please follow the instructions below to let 
# the program know what you want to run and if previously existing files should be used.
# If this is left to the default, all steps will be ran from the start.
#
#   1. Perform quality control on the raw reads.
#      I will run fastqc on the raw reads for each sample and create a folder called QC
#      with the output. 
#      If you did this step already and would like me to skip it, set the run_qc_raw variable
#      below to "No":
run_qc_raw="Yes" # Yes or No
#
#
#   2. Trim the raw reads.
#      I will run trimmomatic on the raw reads to produce high quality alignment-ready reads.
#      If you did this step already and would like to skip it, set the run_trimming variable
#      below to "No"
run_trimming="No" # Yes or No --> This doesn't work - "trimmomatic: command not found" - ask Francesca
#      If you ran trimmomatic yourself, the trimmed reads need to be either in the same folder
#      as the raw reads or in a folder within the project folder called trimmed_reads/, inside 
#      of the corresponding sample subfolder. Look at the directory tree outline if this is unclear.
#      If your reads are along with the raw reads, change the varibable "trimmed_with_raw" from "No" to "Yes"
#      If they have their own folder, leave it at "No".
trimmed_with_raw="Yes"
#
#
#   3. Perform quality control on the trimmed reads.
#      I will also run fastqc on the trimmed reads, for you to check the trimming produced
#      a high enough quality output.
#      If you don't want me to do this, set the run_qc_trimmed variable below to "No"
run_qc_trimmed="Yes" # Yes or No
#
#
#   3. Download the reference transcriptome.
#      I will download the reference transcriptome for Kallisto. If you already have this file
#      locally, uncomment the line below and specify the path to the file.
#transcriptome=<type the path here> 
#      This is the link to the NCBI reference genome I will use. Change it if you'd like to 
#      use a different transcriptome:
transcriptome_link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/405/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_rna.fna.gz"
#
#
#   4. Index the reference transcriptome.
#      I will index the transcriptome for Kallisto to use. 
#      Did you index it already? If so, uncomment below and specify the path to the index file.
#index=<type the path here>
#      Otherwise, comment the line and I will do it for you.
#
#
#   5. Create a translation table connecting transcript names with gene names.
#      Because kallisto works with transcripts, and we are usually interested in genes,
#      we want to add a column to the output of kallisto with gene names. Using the script
#      make_translation.py, we can generate this translation table. 
#      Did you generate it already? If so, uncomment below and specify the path to the json.
#translation_table="/home/projects/dtu_00001/data/TORSK-CELIA/transcriptome/translation_table.json"
#      Otherwise, comment the line and I will do it for you.
### ------------------------------------------------------------------------------------ ###




### --------------------------------- WDIR AND MODULES --------------------------------- ###
## Welcome
echo "Welcome to kallisto!"
echo "This pipeline takes multiple sequencing samples, performs quality checks, trimming and runs the \
pseudoaligner kallisto (https://pachterlab.github.io/kallisto/) on them in a parallel manner."
echo "In this log you will see the steps that were followed along the way."
echo "------"

## Go to the project folder
cd $project_folder
echo "... Project folder is $project_folder"

# # Create output directory (if it doesn't exist already)
# if [ -d $outdir ]
# then 
#     echo "ERROR: The output directory you specified already exists and would be overwritten. \
# Please specify another output directory or remove the existing one."
#     exit 1
# else
#     mkdir $outdir
#     echo "... Created output directory $project_folder$outdir"
# fi
 
# Load all required modules for the job
module load tools
module load kallisto/0.46.0 
module load parallel/20210722
module load anaconda3/4.4.0

module load fastqc/0.11.2
#module load samtools/1.10
module load java/1.8.0
module load trimmomatic/0.38


# module load gcc
# module load intel/perflibs
# module load R/4.1.0
# module load pandoc/2.14.2
### ------------------------------------------------------------------------------------ ###




### ------------------------- QUALITY CONTROL & TRIMMING ------------------------------- ###
# Does the user want to run QC on the raw reads?
echo "------"
if [ $run_qc_raw == "Yes" ] # Yes
then
    echo "... Preparing to run quality control (with fastqc) on the raw reads."
    
    # Here are the raw reads
    R1=$project_folder"raw_reads/S*/S*_1.fq.gz"
    R2=$project_folder"raw_reads/S*/S*_2.fq.gz"

    # Here will go the output of fastqc
    qc_raw_dir=$project_folder"QC/raw_reads/"
    mkdir -p $qc_raw_dir # This might be problematic if the directory already exists

    # Run fastqc (In parallel. Fastqc has the built-in option of parallelizing with -t, so will use that)
    /usr/bin/time -v fastqc -o $qc_raw_dir -t 20 $R1 $R2   

    # Multiqc doesn't work on the queueing system, so the user must run it externally

    echo "... Successfully finished the quality control on the raw reads. I can't run multiqc in the queueing system, so you must run it yourself in the login node. It takes less than a minute."
else # No
    echo "You did not want to run fastqc on the raw reads, so I am skipping this step."
fi



# Does the user want to trim the raw reads?
echo "------"
if [ $run_trimming == "Yes" ] # Yes
then
    # This, at the moment, doesn't work - "trimmomatic: command not found" q
    echo "... Preparing to trim the raw reads (with trimommatic)"

    # Output folder
    trim_dir="trimmed_reads/"
    mkdir -p $trim_dir # Might be problematic if the folder already exists

    # Prepare the list of commands we want to execute in parallel
    commands=()
    samples=()

    cd raw_reads/
    for sample in S*/ # Enter each sample folder within the raw_reads folder
    do
        sample_name=$(echo $sample | cut -d "/" -f 1) # remove the /
        
        # Inputs
        R1="raw_reads/"$sample_name/$sample_name"_1.fq.gz" # Forward raw fastq file
        R2="raw_reads/"$sample_name/$sample_name"_2.fq.gz" # Reverse raw fastq file
        
        # Outputs
        FP=$trim_dir$sample_name/$sample_name"_forward_paired.fq.gz" # Forward paired
        FU=$trim_dir$sample_name/$sample_name"_forward_unpaired.fq.gz" # Forward unpaired
        RP=$trim_dir$sample_name/$sample_name"_reverse_paired.fq.gz" # Reverse paired
        RU=$trim_dir$sample_name/$sample_name"_reverse_unpaired.fq.gz" # Reverse unpaired
        
        cmd="trimmomatic PE -phred33 $R1 $R2 $FP $FU $RP $RU SLIDINGWINDOW:4:15 MINLEN:36"
        commands+=("$cmd")
        samples+=("$sample")
    done

    # Go back to project folder
    cd ..

    # Make sample folders
    cd $trim_dir
    for sample in "${samples[@]}"; do mkdir -p $sample; done
    cd ..

    # Run trimommatic in parallel. Number of parallel runs will equal number of processors
    # specified in the PBS settings (under #PBS -l nodes=1:ppn=20; in this case, 20).
    for command in "${commands[@]}"; do echo $command; done | parallel {}
    
    echo "... Successfully finished trimming the reads."
else # No
    echo "You did not want to run trimmomatic, so I am skipping this step and will use the existing files."
fi




# Does the user want to run QC on the trimmed reads?
echo "------"
if [ $run_qc_trimmed == "Yes" ] # Yes
then
    echo "... Preparing to run quality control (with fastqc) on the trimmed reads."
    
    if [ $trimmed_with_raw == "Yes" ] # Trimmed reads are with raw reads
    then
        # Here are the trimmed reads
        R1=$project_folder"raw_reads/S*/*_forward_paired.fq.gz"
        R2=$project_folder"raw_reads/S*/*_reverse_paired.fq.gz"

        # Here will go the output of fastqc
        qc_trimmed_dir=$project_folder"QC/trimmed_reads/"
        mkdir -p $qc_trimmed_dir # This might be problematic if the directory already exists

        # Run fastqc (In parallel. Fastqc has the built-in option of parallelizing with -t, so will use that)
        /usr/bin/time -v fastqc -o $qc_trimmed_dir -t 20 $R1 $R2 
        
    else # Trimmed reads are in their own folder
        # Here are the trimmed reads
        R1=$project_folder$trim_dir"/S*/*_forward_paired.fq.gz"
        R2=$project_folder$trim_dir"/S*/*_reverse_paired.fq.gz"

        # Here will go the output of fastqc
        qc_trimmed_dir=$project_folder"QC/trimmed_reads/"
        mkdir -p $qc_trimmed_dir # This might be problematic if the directory already exists

        # Run fastqc (In parallel. Fastqc has the built-in option of parallelizing with -t, so will use that)
        /usr/bin/time -v fastqc -o $qc_trimmed_dir -t 20 $R1 $R2
        
    fi

    # Multiqc doesn't work on the queueing system, so the user must run it externally

    echo "... Successfully finished the quality control on the trimmed reads. I can't run multiqc in the queueing system, so you must run it yourself in the login node. It takes less than a minute."

else # No
    echo "You did not want to run fastqc on the trimmed reads, so I am skipping this step."
fi
### ------------------------------------------------------------------------------------ ###




### ----------------------- REFERENCE TRANSCRIPTOME + KALLISTO INDEX -------------------- ###
# If the transcriptome path is not specified, download it and set filename variable
if [ -z "$transcriptome" ]  
then 
    echo "... A reference transcriptome file was not provided, so it will be downloaded \
from the internet at $transcriptome_link"
    wget $transcriptome_link --directory-prefix $project_folder 
    transcriptome=$project_folder$(echo $transcriptome_link | rev | cut -d '/' -f1 | rev)
    echo "... The transcriptome was successfully downloaded! It can be found at $transcriptome"
else
    echo "... The reference transcriptome file was provided at $transcriptome."
fi


# If the index file is not specified, generate it and set filename variable
if [ -z "$index" ]
then
    echo "... Indexing the reference transcriptome."
    index=$project_folder"index"
    /usr/bin/time -v kallisto index -i $index $transcriptome
    echo "... Successfully indexed! Find index at $index."
else
    echo "... Indexed transcriptome was provided at $index. "
fi


# If the translation table is not specified, generate it and set filename variable
if [ -z "$translation_table" ]
then 
    echo "... Making transcript-gene name translation table."
    translation_table=$project_folder"translation_table.json"
    /usr/bin/time -v python3 $make_translation $transcriptome $translation_table
    echo "... Successfully generated the table. Find it at $translation_table"
else
    echo "... Transcript-gene name translation table was provided at $translation_table"
fi
### ------------------------------------------------------------------------------------ ###




### ---------------------------------- KALLISTO QUANT ----------------------------------- ###
echo "------"
echo "... Preparing to run kallisto quantification algorithm."
# Prepare the list of commands we want to execute in parallel
commands=()
samples=()

# Make output folder
mkdir quantification_kallisto/

# To trimmed reads
if [ $trimmed_with_raw == "Yes" ] # Trimmed reads are with raw reads
then
    cd raw_reads/
    for sample in S*/ # Enter each sample folder within the wdir
    do
        F1="raw_reads/"$sample*_forward_paired.fq.gz # Trimmed forward paired
        F2="raw_reads/"$sample*_reverse_paired.fq.gz # Trimmed reverse paired
        out="quantification_kallisto/"$sample 
        cmd="/usr/bin/time -v kallisto quant -i $index $F1 $F2 -o $out"
        commands+=("$cmd")
        samples+=("$sample")
    done
else # Trimmed reads are in their own folder
    cd $trim_dir
    for sample in S*/ # Enter each sample folder within the wdir
    do
        F1=$trim_dir$sample*_forward_paired.fq.gz # Trimmed forward paired
        F2=$trim_dir$sample*_reverse_paired.fq.gz # Trimmed reverse paired
        out="quantification_kallisto/"$sample 
        cmd="/usr/bin/time -v kallisto quant -i $index $F1 $F2 -o $out"
        commands+=("$cmd")
        samples+=("$sample")
    done
fi



# Back to project folder
cd ..

# Run kallisto in parallel. Number of parallel runs will equal number of processors
# specified in the PBS settings (under #PBS -l nodes=1:ppn=20; in this case, 20).
for element in "${commands[@]}"; do echo $element; done | parallel {}

echo "... Successfully ran kallisto on the following samples:"
for sample in "${samples[@]}"; do echo $sample; done
echo "Results can be found in the directories named after the samples."
### ------------------------------------------------------------------------------------ ###




### ---------------------- AUGMENT KALLISTO OUTPUT WITH GENE NAMES --------------------- ### 
echo "------"
echo "... Preparing to augment kallisto output with the corresponding gene names."
cd quantification_kallisto/

for sample in S*/
do
    abundance_file=$sample"abundance.tsv"
    output_file=$sample"abundance_aug.tsv"
    /usr/bin/time -v python3 $add_gene_name $abundance_file $output_file $translation_table
done

echo "... Successfully finished the augmentation."
### ------------------------------------------------------------------------------------ ###




### -------------------------------- MERGED GENE COUNTS -------------------------------- ###
# For the DESeq2 analysis, we need the results of kallisto in the format:
#       sample1 sample2 ... sampleN
# gene1    25      130         34
# gene2     4        0          7  
# ...
# geneM   124        3        101

# In addition, because kallisto works with transcripts and not genes, 
# multiple transcripts may map to the same gene and rows must be collapsed. 
# We decide to do so by adding together the number of reads that map to the same
# gene (but different isoforms), because our analysis shows it is the most
# accurate way of doing it - See Correlation-DealingWithIsoforms.ipynb

merged_file=$project_folder"merged_samples.tsv"

/usr/bin/time -v python3 $merge_samples $project_folder"quantification_kallisto/" $merged_file
### ------------------------------------------------------------------------------------ ###




### ------------------------------------------------------------------------------------ ###
echo "------"
echo "Congratulations! Your pipeline has completed."
echo "You can find your merged results in $merged_file. Please check the error log kallisto.err to see if there were any issues."
echo "You can use the RMarkdown notebook in $deseq2notebook to perform the differential expression analysis on the output of kallisto. The suggested way of doing this is exporting both the merged_samples.tsv file and the RMarkdown to your local computer and run the notebook interactively in RStudio."
echo "Good luck with it!"
echo "------"
### ------------------------------------------------------------------------------------ ###