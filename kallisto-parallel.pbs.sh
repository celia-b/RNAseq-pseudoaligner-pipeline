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
## SAMPLES
# Each sample must have a folder called S+number (e.g. S1, S2, S100...)
# Inside of these folders, the raw/trimmed fastq files must be located. 
#
# Please state the path containing those sample folders:
workdir="/home/projects/dtu_00001/data/TORSK/"


## OUTPUT
# The output of kallisto will be written in the wdir, inside of a 
# user-defined folder. Each sample will have its own folder, with the
# same names as the fastq folders (e.g. S1, S2, S100...)
#
# Please specify the name of the folder containing kallisto output for all samples:
outdir="kallisto_test/"


## PREREQUIREMENTS
# Kallisto has a series of requirements that must either be ran beforehand and pointed
# to for this script to use, or can be executed directly by this script.
#
# These requirements are:
#   1. Download the reference transcriptome.
#      a. Did you download it already? If so, uncomment below and specify the path to it: 
#transcriptome="/home/projects/dtu_00001/data/TORSK-CELIA/transcriptome/Gadus_morhua_NCBI.fa.gz" 
#      b. Otherwise, comment the line above and we will do it for you.
#         Here you can update the link to the NCBI reference genome you want to use:
transcriptome_link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/405/GCF_902167405.1_gadMor3.0/GCF_902167405.1_gadMor3.0_rna.fna.gz"
#
#   2. Index the reference transcriptome.
#      a. Did you index it already? If so, uncomment below and specify the path to the index.
#index="/home/projects/dtu_00001/data/TORSK-CELIA/transcriptome/ncbi_kallisto_index"
#      b. Otherwise, comment the line and we will do it for you.
#
# 3. Create a translation table connecting transcript names with gene names.
#    Because kallisto works with transcripts, and we are usually interested in genes,
#    we want to add a column to the output of kallisto with gene names. Using the script
#    make_translation.py, we can generate this translation table. 
#    a. Did you generate it already? If so, uncomment below and specify the path to the json.
#translation_table="/home/projects/dtu_00001/data/TORSK-CELIA/transcriptome/translation_table.json"
#    b. Otherwise, comment the line and we will do it for you.


# Temporary specification of where make_translation.py and add_gene_name.py are located.
# Should instead be in the same folder as this script.
make_translation="/home/projects/dtu_00001/data/TORSK-CELIA/kallisto_pipeline/make_translation.py"
add_gene_name="/home/projects/dtu_00001/data/TORSK-CELIA/kallisto_pipeline/add_gene_name.py"
### ------------------------------------------------------------------------------------ ###




### --------------------------------- WDIR AND MODULES --------------------------------- ###
# Welcome
echo "Welcome to kallisto!"
echo "This pipeline takes multiple sequencing samples and runs the \
pseudoaligner kallisto (https://pachterlab.github.io/kallisto/) on them in a parallel manner."
echo "In this log you will see the steps that were followed along the way."
echo "------"

# Go to working directory
cd $workdir
echo "... Working directory is $workdir"

# Create output directory (if it doesn't exist already)
if [ -d $outdir ]
then 
    echo "ERROR: The output directory you specified already exists and would be overwritten. \
Please specify another output directory or remove the existing one."
    exit 1
else
    mkdir $outdir
    echo "... Created output directory $workdir$outdir"
fi
 
# Load all required modules for the job
module load tools
module load kallisto/0.46.0 
module load parallel/20210722
module load anaconda3/4.4.0
### ------------------------------------------------------------------------------------ ###




### ----------------------- REFERENCE TRANSCRIPTOME + KALLISTO INDEX -------------------- ###
# If the transcriptome path is not specified, download it and set filename variable
if [ -z "$transcriptome" ]  
then 
    echo "... A reference transcriptome file was not provided, so it will be downloaded \
from the internet at $transcriptome_link"
    wget $transcriptome_link --directory-prefix $workdir$outdir 
    transcriptome=$workdir$outdir$(echo $transcriptome_link | rev | cut -d '/' -f1 | rev)
    echo "... The transcriptome was successfully downloaded! It can be found at $transcriptome"
else
    echo "... The reference transcriptome file was provided at $transcriptome."
fi


# If the index file is not specified, generate it and set filename variable
if [ -z "$index" ]
then
    echo "... Indexing the reference transcriptome."
    index=$workdir$outdir"index"
    /usr/bin/time -v kallisto index -i $index $transcriptome
    echo "... Successfully indexed! Find index at $index."
else
    echo "... Indexed transcriptome was provided at $index. "
fi


# If the translation table is not specified, generate it and set filename variable
if [ -z "$translation_table" ]
then 
    echo "... Making transcript-gene name translation table."
    translation_table=$workdir$outdir"translation_table.json"
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

for dir in S*/ # Enter each sample folder within the wdir
do
    F1=$dir*_forward_paired.fq.gz # Forward fastq file
    F2=$dir*_reverse_paired.fq.gz # Reverse fastq file
    out=$outdir$dir
    cmd="/usr/bin/time -v kallisto quant -i $index $F1 $F2 -o $out"
    commands+=("$cmd")
    samples+=("$dir")
done

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
cd $workdir$outdir

for dir in S*/
do
    abundance_file=$dir"abundance.tsv"
    output_file=$dir"abundance_aug.tsv"
    /usr/bin/time -v python3 $add_gene_name $abundance_file $output_file $translation_table
done

echo "... Successfully finished the augmentation."
### ------------------------------------------------------------------------------------ ###




### ------------------------------------------------------------------------------------ ###
echo "------"
echo "Congratulations! Your pipeline has completed."
echo "Please check the error log kallisto.err to see if there were any issues."
echo "Good luck with your analysis!"
### ------------------------------------------------------------------------------------ ###