#!/usr/bin/env python
# coding: utf-8

# In[92]:


import pandas as pd
import subprocess
import os

full_tab = pd.read_table('Data/Sup_tables/sup_table_3_single_chain.tsv',sep='\t')

sifts_tab = pd.read_table("path_to_SIFTS_data",sep="\t", skiprows=1)
sifts_tab.columns = ['pdb', 'chain', 'sp_primary', 'res_beg', 'res_end', 'pdb_beg', 'pdb_end', 'sp_beg', 'sp_end']


ecod_tab = pd.read_csv("path_to_ECOD_data", skiprows=4, header=0, sep='\t')

# keep only ecod families that are within the af2 non redundant set, "full_tab"
ecod_tab_of_relevant_families = ecod_tab[ecod_tab.f_id.isin(full_tab.f_id.unique())]
# get all the pdbs of these families based on ecod
pdbs_of_families = ecod_tab_of_relevant_families.pdb.unique()

# keep only sifts entries for the pdbs matching ecod families from af2 non redundant
filtered_sifts = sifts_tab[sifts_tab.pdb.isin(pdbs_of_families)]

# filt_sifts_grpby_pdb_and_chain = filtered_sifts.groupby(['pdb','chain']).size()
# multi_row_pdbs = filt_sifts_grpby_pdb_and_chain[filt_sifts_grpby_pdb_and_chain > 1]

# filter duplicate rows - any gapped (unresolved) segment in a pdb induces a new line in sifts 
filtered_no_dup = filtered_sifts.drop_duplicates(subset=['pdb','chain','sp_primary'])

# group the filtered sifts tab by pdb and chain, compute the size of each group
grp_sizs = filtered_no_dup.groupby(['pdb','chain']).size()
# keep only groups of size >=2, meaning more than 1 uniprot for a chain; 
# merge with the group series with the filtered sifts data
multi_uniprot_tab = grp_sizs[grp_sizs >= 2].reset_index().drop(0,axis=1).merge(filtered_no_dup, how='left',on=['pdb','chain'])

# # merge with the ecod data 
# # Note that NaNs could occure for cases where a chain doesn't have a mapped ecod family, e.g. a peptide chain. 
# multi_uniprot_with_families = multi_uniprot_tab.merge(ecod_tab_of_relevant_families[['f_id','pdb','chain']], how='left',on=['pdb','chain'])

# check that indeed all families are also in the af2 non redun tab, after dropping NaNs
sum(multi_uniprot_with_families.dropna(subset=['f_id']).f_id.isin(full_tab.f_id.unique())) == len(multi_uniprot_with_families.dropna(subset=['f_id']))

multi_uniprot_matched_families = multi_uniprot_tab.merge(ecod_tab_of_relevant_families[['f_id','pdb','chain','unp_acc']].rename(columns={'unp_acc':'sp_primary'}))
multi_uniprot_matched_families.merge(full_tab, on=['f_id']).to_csv('multi_uniprot_matched_families.tsv',sep="\t", index=False)

