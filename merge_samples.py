#!usr/bin/env python3
import sys
import os
import re

import pandas as pd
from functools import reduce

import gc


## -------- User input -------- ##
## Directory containing all sample folders with kallisto results & output file
try:
    kallisto_dir = sys.argv[1] 
    samples = [s for s in os.listdir(kallisto_dir) if s.startswith('S')]
    kallisto_files = [kallisto_dir + s + '/' + f for s in samples for f in os.listdir(kallisto_dir + s) if f.startswith('abundance_aug')]
    output = sys.argv[2]
except:
    print('ERROR: This program requires two arguments: (1) directory containing all sample folders with kallisto results (2) desired output file.')
    sys.exit(1)
## ---------------------------- ##




## --- Merge kallisto files --- ##
print('Loading files...')
lst_dfs = []
for file in kallisto_files:
    try:
        df = pd.read_csv(file, sep='\t', header=0, usecols=['target_id', 'gene_name', 'est_counts'])
        df = df[['target_id', 'gene_name', 'est_counts']] # Sort columns
        sample = re.search('\/(S.*)\/abundance_aug.tsv$', file).group(1) # Get sample name
        df.rename(columns = {'est_counts': sample}, inplace=True) # Name counts column after sample
        lst_dfs.append(df)
        
        # Clear up mem
        del df
        gc.collect()
    
    except pd.errors.EmptyDataError:
        print(file + ' is empty.')
        continue

print('Merging files...')
kallisto_all = reduce(lambda left, right: pd.merge(left, right, how='inner', on=['target_id', 'gene_name']), lst_dfs)

# Clear up mem
del lst_dfs
gc.collect()
## ---------------------------- ##




## ---- Aggregate isoforms ---- ##
# Kallisto works with transcripts, when we usually study genes. 
# For this reason, we decide to aggregate the counts from rows (transcripts)
# that correspond to the same gene.

# NOTE: If we wanted to study isoforms separately, this step would
#       have to be ommitted.

print('Aggregating isoforms...')
kallisto_agg = kallisto_all.groupby('gene_name').agg('sum')
## ---------------------------- ##




## ---- Write merged file ----- ##
print('Writing to file...')
kallisto_agg.to_csv(output, sep='\t')
## ---------------------------- ##

