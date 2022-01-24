#!usr/bin/env python3
import gzip
import re
import json
import sys

## ----- make_translation.py ----- ##
# This program makes a transcript-to-gene name translation dictionary 
# from transcriptome fasta file.

# In transcriptome fasta files, the format of the headers is:
# >transcript_name description (gene_name), type_of_sequence
# For example:
# >XM_030337033.1 PREDICTED: Gadus morhua family with sequence similarity 43 member B (fam43b), mRNA
# With regular expressions, we can extract the transcript_name and gene_name 
# and make a translation dictionary with the equivalences
## ------------------------------ ##


## ------ User input ------ ##
# abundance.tsv file, output file and transcript-gene name translation table.
try:
    transcriptome_file = sys.argv[1]
    output = sys.argv[2]
except:
    print('ERROR: You must specify the full path of your reference transcriptome file (gzipped) \
as first argument and of the output file as second.')
    sys.exit(1)
## ------------------------ ##


## --------- Main ---------- ##
translation = {}
with gzip.open(transcriptome_file, "rb") as transcriptome:
    for line in transcriptome:
        if line.startswith(b'>'):
            # Find transcript and gene name pairs
            transcript = re.search(b'^>(\w*\.*\w*)\s', line)
            gene = re.findall(b'\((\w+-*\w*)\),', line)
            
            # Add transcript:gene as key:value pair in dict
            translation[transcript.group(1).decode('ascii')] = gene[-1].decode('ascii')

            # Debug: Check if something other than the transcript and gene names are being captured by the regex
            try:
                print(transcript.group(2))
                print(gene.group(2))
            except Exception:
                pass
## ------------------------- ##



## -------- Output --------- ##
# Write json file from translation dict.
# This output file can later be used to augment the Kallisto abundances.tsv files.
# To do that, use python program called 'add_gene_name.py'

with open(output, 'w') as json_file:
    json.dump(translation, json_file)
## ------------------------- ##
