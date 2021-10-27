#!/usr/bin/python3

# script that gets models (af2 or others) and computes secondary strutrue for them using DSSP method

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import pandas as pd
import glob
import re

# set dssp bin path
path_to_dssp = "local_path_to_dssp"

# helix and strand length definition
helix_length = 3
strand_length = 3

# models to use, replace "my_model" with relevant prefix
pdb_for_dssp = glob.glob('*my_model*')[0]


p = PDBParser()
structure = p.get_structure("id", pdb_for_dssp)
model = structure[0]
dssp = DSSP(model, pdb_for_dssp, dssp=path_to_dssp)

dssp_df = pd.DataFrame.from_dict(dssp)
dssp_df.rename(columns={2:'ss2'},inplace=True)

dssp_str = dssp_df['ss2'].str.cat(sep='')

# replace 8 letter code with 3 letter code
def replace_8letter_code(ss2):
    ss2 = re.sub('[GI]','H',ss2)
    ss2 = re.sub('[B]','E',ss2)
    ss2 = re.sub('[ST]','C',ss2)
    ss2 = re.sub('[ ]','C',ss2)
    ss2 = re.sub('-','C',ss2)
    return ss2

# define ss2 groups based on helix and strand definitions
def define_ss2_groups(ss2):
    if ss2 != 'ss2_missing':

        helix_param = 'H'*helix_length
        strand_param = 'E'*strand_length
        is_helix = helix_param in ss2
        is_strand = strand_param in ss2
       
        if is_helix and is_strand:
            is_combined = True
        else:
            is_combined = False
        
    else:
        is_helix = is_strand = is_combined = None
        
    return (is_helix, is_strand, is_combined)

# run on the input structure
dssp_str = replace_8letter_code(dssp_str)
dssp_groups = define_ss2_groups(dssp_str)

# save outputs
dssp_results_id = re.sub('for_dssp','dssp_results',pdb_for_dssp) # chane the input file names according to relevant convention in directory
dssp_results_id = re.sub('.pdb','.tsv',dssp_results_id)

merged_output = (*dssp_groups, dssp_str)
dssp_row = pd.DataFrame(pd.Series(merged_output)).transpose()
dssp_row.rename(columns={0:'pep_is_helix',1:'pep_is_strand',2:'pep_is_combined',3:'pep_ss2'}, inplace = True)

dssp_row.to_csv(dssp_results_id, sep='\t',index=False)

