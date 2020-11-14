#!/usr/bin/env python3

import io
import os
import pandas as pd
import numpy as np
import argparse as ap


parser = ap.ArgumentParser(description="script to split sample information in vcf columns")
parser.add_argument('a', help='vcf file to split')
parser.add_argument('b', help='output_vcf_file')
args = parser.parse_args()

fileout = open(args.b, 'w')

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


vcf_input = read_vcf(args.a)

# extract longest column name for list
field_length = vcf_input.FORMAT.astype(str).map(len)
longest_col = vcf_input.loc[field_length.argmax(), 'FORMAT']
longest_col_list = longest_col.split(':')

# rename columns to include sample information (tumour or Normal) 
normal_cols = [s + ' Normal' for s in longest_col_list]
tumour_cols = [s + ' Tumour' for s in longest_col_list]

# split column informaiton into new columns for both tumour and normal
vcf_input = vcf_input.join(vcf_input.iloc[:, 9].str.split(':', expand=True).add_prefix('Tumour sample'))
vcf_input = vcf_input.join(vcf_input.iloc[:, 10].str.split(':', expand=True).add_prefix('Normal sample'))

# get column names from data frame and get first 10 columns
vcf_input_cols = vcf_input.columns.tolist()
vcf_input_cols = vcf_input_cols[0:11]

# create list of all columns with new headers 
vcf_input_cols = vcf_input_cols + tumour_cols + normal_cols

# replace col names with new ones from list
vcf_input.columns = vcf_input_cols

vcf_input.to_csv(fileout, header=True, index=False, sep='\t')

