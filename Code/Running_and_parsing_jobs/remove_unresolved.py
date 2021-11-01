#!/usr/bin/env python
# coding: utf-8

# This script aligns models to the native and removes regions that are unresolved in the native structures.
# This is required for calculating measure with FlexPepDock, which cannot handle model and native structures
# with differing length. Input is directory with the models and a native file.

# In[28]:


from biopandas.pdb import PandasPdb
from Bio import pairwise2
import pandas as pd
import glob
import os
import sys


# In[22]:


def get_chain_seq(pdb, chain):
    sequence = pdb.amino3to1()
    sequence_list = list(sequence.loc[sequence['chain_id'] == chain, 'residue_name'])
    seq = ''.join(sequence_list)

    return(seq)


# In[23]:


def remove_non_resolved_from_a_chain(chain_native, chain_model, native, model):
    chain_seq_native = get_chain_seq(native, chain_native)
    chain_seq_model = get_chain_seq(model, chain_model)
    aln = pairwise2.align.globalxs(chain_seq_native, chain_seq_model, -3, -1) # align with high penalty for gaps, to keep especially the peptide together
    aln_native = aln[0][0]
    aln_model = aln[0][1]

    res_numbers = []
    for l in range(0, len(aln_native)):
        if(aln_native[l] != '-'):
            res_numbers.append(l+1)

    model_new_df = model.df['ATOM'][(model.df['ATOM']['residue_number'].isin(res_numbers)) & (model.df['ATOM']['chain_id'] == chain_model)]
    return(model_new_df)


# In[24]:


def renumber_peptide_chain(model):
    first_residue_in_chain = min(model.df['ATOM'].query('chain_id == "B"')['residue_number'])
    model.df['ATOM'].loc[(model.df['ATOM'].chain_id == 'B'),'residue_number'] = model.df['ATOM'].loc[(model.df['ATOM'].chain_id == 'B'),'residue_number'] - first_residue_in_chain + 1
    return(model)


# In[25]:


def truncate_a_model(chain_rec, chain_pep, native, model_file):
    model=PandasPdb().read_pdb(model_file)
    model = renumber_peptide_chain(model)

    rec_df = remove_non_resolved_from_a_chain(chain_rec, 'A', native, model)
    pep_df = remove_non_resolved_from_a_chain(chain_pep, 'B', native, model)
    new_atoms = pd.concat([rec_df, pep_df], ignore_index=True)
    model_new = model
    model_new.df['ATOM'] = new_atoms

    model_new_file = model_file.replace('.pdb', '_truncated.pdb')
    print(model_new_file)
    model_new.to_pdb(path=model_new_file, records=['ATOM'], gz=False, append_newline=False)


# In[26]:


dir = sys.argv[1]
native_file = sys.argv[2]
complex = os.path.basename(os.path.dirname(dir)).split('_')

chain_rec = 'A'
chain_pep = 'B'

all_models = glob.glob(os.path.join(dir, '*superimposed.pdb'))

native = PandasPdb().read_pdb(native_file)

for model_file in all_models:
    truncate_a_model(chain_rec, chain_pep, native, model_file)


