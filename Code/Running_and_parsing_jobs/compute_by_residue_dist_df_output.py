#!/usr/bin/python3

# script to compute the rmsd values residue-by-residue compared to the native complex. 
# this script also computes the by residue LDDT, and outputs rmsd table, LDDT table and a comparison table between both. 
# script assums running in the directory with the models and native

from biopandas.pdb import PandasPdb
import pandas as pd
import os 
import os.path
import pandas as pd 
import Bio
from Bio import SeqIO
import pickle 
import re
import sys
from Bio import pairwise2
import copy 
import glob
import numpy as np

# provide the id of native_pdb 
native_pdb_id = sys.argv[1]
print(native_pdb_id)

# function to process and prepare the native into a dictionary of chains ids, sequences (rec and pep), lengths etc. 

def prepare_native(native_pdb):
    native_dict = {}
    npdb = PandasPdb()
    npdb.read_pdb(native_pdb)
    native_df=npdb.df['ATOM']
    native_dict['native_df'] = native_df.copy(deep=True)
    
    native_chains=native_df.chain_id.unique()
    first_chain_len=len(native_df.query('chain_id == @native_chains[0]'))
    second_chain_len=len(native_df.query('chain_id == @native_chains[1]'))

    if first_chain_len > second_chain_len:
        rec_chain_native = native_chains[0]
        pep_chain_native = native_chains[1]
    else:
        rec_chain_native = native_chains[1]
        pep_chain_native = native_chains[0]
    
    native_dict['rec_chain_native'] = rec_chain_native
    native_dict['pep_chain_native'] = pep_chain_native
    
    rec_seq_native=npdb.amino3to1().query('chain_id == @rec_chain_native')
    pep_seq_native=npdb.amino3to1().query('chain_id == @pep_chain_native')
    rec_seq_native="".join(list(rec_seq_native.residue_name))
    pep_seq_native="".join(list(pep_seq_native.residue_name))
  
    
    native_pep_df = native_df.query('chain_id == @pep_chain_native')
    native_first_pos = native_pep_df.residue_number.min()
    native_dict['native_pep_df'] = native_pep_df.copy(deep=True)
    native_dict['native_first_pos'] = native_first_pos
    
    native_dict['rec_seq_native'] = rec_seq_native
    native_dict['pep_seq_native'] = pep_seq_native
    
    return native_dict 

# function to process and prepare each model into a dictionary of chains ids, sequences (rec and pep), lengths etc. 

def prepare_models_and_compute_distance(model_pdb, dict_of_native):
        
    mpdb = PandasPdb()
    mpdb.read_pdb(model_pdb)
    model_df=mpdb.df['ATOM']
    
    rec_seq_model=mpdb.amino3to1().query('chain_id == "A"')
    pep_seq_model=mpdb.amino3to1().query('chain_id == "B"')
    rec_seq_model="".join(list(rec_seq_model.residue_name))
    pep_seq_model="".join(list(pep_seq_model.residue_name))

    pep_alignment = pairwise2.align.globalxs(dict_of_native['pep_seq_native'], pep_seq_model, -3, -1)
    
    model_pep_df = model_df.query('chain_id == "B"')
    model_first_pos = model_pep_df.residue_number.min()
    
    native_offset = 0
    model_offset = 0
    rms_list = []
    pos_list = []
    for i in range(len(pep_alignment[0][1])): # the string of the model peptide fasta 
        if pep_alignment[0][0][i]=="-":
#            print("skipping")
            model_offset+=1
            rms_list.append(-1) # not computed
        elif pep_alignment[0][1][i]=="-":
#            print("skipping")
            native_offset+=1
            rms_list.append(-1)
        else:
            rms_val, pos_pair = by_residue_pair_dist(native_offset, model_offset, dict_of_native['native_pep_df'],
                                 dict_of_native['native_first_pos'], model_pep_df, model_first_pos)
            rms_list.append(rms_val)
            pos_list.append(pos_pair)



#            rms_list.append(by_residue_pair_dist(native_offset, model_offset, dict_of_native['native_pep_df'],
#                                 dict_of_native['native_first_pos'], model_pep_df, model_first_pos))
            model_offset+=1
            native_offset+=1
    pos_name = str(native_pdb_id[0:4])+"_position_list.txt"
    pd.Series(pos_list).to_csv(pos_name)
    lddt_series = model_df.query('chain_id == "B"').drop_duplicates(subset='residue_number').b_factor
    return pd.Series(rms_list), lddt_series, pep_seq_model

# helper function to compute the offset of residue numbering between model and native structures
def by_residue_pair_dist(native_offset, model_offset, native_pep_df, native_first_pos, model_pep_df, model_first_pos):
    
    
    native_i_pos = native_first_pos + native_offset
    model_i_pos = model_first_pos + model_offset
    filt_native_pep_df, filt_model_pep_df = intersect_res_dfs(
        native_pep_df.query('residue_number == @native_i_pos'),
        model_pep_df.query('residue_number == @model_i_pos')
    )
    i_rms = PandasPdb.rmsd(
    filt_native_pep_df,
    filt_model_pep_df
    )
    one_based_model_offset = model_offset+1
    return i_rms, (native_i_pos, model_i_pos, one_based_model_offset)


# only keep the atoms that overlap in each compared residue (e.g. model has OXT and native doesnt)
def intersect_res_dfs(df1, df2):
    filt_df1 = df1[df1.atom_name.isin(df2.atom_name)]
    filt_df2 = df2[df2.atom_name.isin(df1.atom_name)]
    
    filt_df1.sort_values('atom_name',inplace=True)
    filt_df2.sort_values('atom_name',inplace=True)
    return filt_df1, filt_df2


# make the output table of rms and lddt
def make_rms_lddt_tab(columns_lst, lddt_series, rms_series, i_model):
    col_series = pd.Series(columns_lst)
    col_series = col_series.reset_index()
    col_series.rename(columns={0:"seq"}, inplace=True)
    lddt_series = lddt_series.reset_index()[i_model]

    rms_series = rms_series.reset_index()[i_model]
    lddt_series.name = "lddt"
    rms_series.name = "rms"


    full=pd.concat([col_series,rms_series,lddt_series],axis=1)
    full = full[['seq','rms','lddt']]
    full['model'] = i_model
    return full

# Main - run the whole thing
native_dict = prepare_native(native_pdb_id)
dict_copy = copy.deepcopy(native_dict)

series_lst =[]
ld_vs_rms_lst = []
lddt_lst = []

model_lst = glob.glob("*my_model*.pdb") # replace the prefix with relevant model name

# loop over all the relevant models
for i in range(len(model_lst)):
    
    first_model_occur=model_lst[i].find("model_")
    row_name=model_lst[i][first_model_occur:(first_model_occur+7)] #+7 gets the "model_5" part of the output name, which is what network params were used (in our naming setup)
    
    dict_copy = copy.deepcopy(native_dict)
    i_series, lddt_series, pep_seq_model =prepare_models_and_compute_distance(model_lst[i], dict_copy)
    columns_lst = [pep_seq_model[x] for x in range(len(pep_seq_model))]

    i_series.name=row_name
    series_lst.append(i_series)
    
    
    lddt_series.name=row_name
    lddt_lst.append(lddt_series)
    
    ld_vs_rms_lst.append(make_rms_lddt_tab(columns_lst,lddt_series,i_series,row_name))

    del dict_copy

# save outputs as csv
by_resi_df = pd.concat(series_lst,axis=1).transpose()
by_resi_df.columns = columns_lst
by_resi_df.to_csv(str(native_pdb_id)[0:4]+'_by_residue_rms.tsv',sep="\t")

lddt_df = pd.concat(lddt_lst, axis=1).transpose()
lddt_df.columns = columns_lst
lddt_df.to_csv(str(native_pdb_id)[0:4]+'_lddt_by_residue.tsv',sep="\t")

rms_lddt_tab = pd.concat(ld_vs_rms_lst)
rms_lddt_tab['pdb'] = str(native_pdb_id)[0:4]
rms_lddt_tab.to_csv(str(native_pdb_id)[0:4]+'_rms_vs_lddt_comparison_by_residue.tsv',sep="\t", index=False)

