#!usr/bin/env python3
import gzip
import re
import json
import sys

## ------ add_gene_name.py ------ ##
# This program uses the output from 'make_translation.py' (a json file containing 
# the transcript-gene name translation dict) to add a column to the Kallisto
# 'abundance.tsv' output file with the gene names corresponding to each transcript id.
## ------------------------------ ##


## ------ User input ------ ##
# abundance.tsv file, output file and transcript-gene name translation table.
try:
    input = sys.argv[1]
    output = sys.argv[2]
    translation_json = sys.argv[3]
except:
    print('ERROR: You must specify the full path of your abundance.tsv file as first argument, \
of your output file as second and of the translation table generated by make_translation.py as third.')
    sys.exit(1)
## ------------------------ ##


## ----- Load translation dict ------- ##
# Read output from 'make_translation.py'
with open(translation_json) as json_file:
    translation_table = json.load(json_file)
## ----------------------------------- ##


## --------- Main + Output --------- ##
# Read 'abundance.tsv' file and find transcript-gene correspondances.
# Write to output file 'abundance_aug.tsv'.

with open(input, 'r') as abundance_file:
    with open(output, 'w') as augmented_file:
        for line in abundance_file:
            line = line.strip()
            # If it is the header line, add column name
            if line.startswith('target_id'):
                line = line + '\tgene_name'
                print(line, file = augmented_file)
            
            # Else:
            else: 
                # Find transcript name
                transcript = re.search('^(\w*\.*\w*)\s', line).group(1)
                
                # Add it to the line, if it is found in the translation dict
                try:
                    gene = translation_table[transcript]
                    line = line + '\t' + gene
                    print(line, file = augmented_file)
                
                # Otherwise just print the line
                except KeyError as err:
                    print(line, file = augmented_file)
                    # Debug
                    print('A transcript name could not be found in the translation dict')
## ------------------------------ ##
