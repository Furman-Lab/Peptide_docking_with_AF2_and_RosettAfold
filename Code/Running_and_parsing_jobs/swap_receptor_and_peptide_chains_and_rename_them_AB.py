#!/usr/bin/env python
# coding: utf-8

from biopandas.pdb import PandasPdb
import pandas as pd
from pandas import Series
import os 
import sys

# pdb to process
pdb_to_read = sys.argv[1]

ppdb = PandasPdb()
ppdb.read_pdb(pdb_to_read)

# keep only ATOM and HETATM rows
df = pd.concat([ppdb.df['ATOM'],ppdb.df['HETATM']])

def order_and_rename(df, first_chain, second_chain):

    first_df = df.query("chain_id == @first_chain")
    first_df.chain_id = "A"
    second_df = df.query("chain_id == @second_chain")
    second_df.chain_id = "B"

    merged_df = pd.concat([first_df, second_df])
    merged_df.line_idx = list(range(len(merged_df)))
    return merged_df

new_pdb = PandasPdb()

chains = list(set(df.chain_id))

# if first chain is peptide, reorder and rename; if first chain is receptor - rename only
if len(df.query('chain_id == @chains[0]')) < len(df.query('chain_id == @chains[1]')):
    new_pdb.df['ATOM']= order_and_rename(df, chains[1], chains[0])
else:
    new_pdb.df['ATOM']= order_and_rename(df, chains[0], chains[1])


# save out pdb
new_pdb.to_pdb(pdb_to_read[0:4]+'_reordered_renamed_native.pdb')

