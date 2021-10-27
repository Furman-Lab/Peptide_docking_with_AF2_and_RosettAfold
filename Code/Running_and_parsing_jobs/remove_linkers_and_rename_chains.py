#!/usr/bin/env python
# coding: utf-8

## script to remove the linkers from af2 outputs that were modeled with a 30GLY linker. also renames the chains to be A for receptor and B for peptide

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

# pdb to process
pdb_to_read = sys.argv[1]

# table with fastas of the modeled pdbs
fasta_tab = sys.argv[2]

# a3m file of msa from af2 output
a3m_file = sys.argv[3]

# pdb id of the model
i_pdb_id = str(sys.argv[4]).lower()

# chains of rec and peptide
chains = str(sys.argv[5])

# is it a polyA run?
polyA = str(sys.argv[6])


record_dict = SeqIO.to_dict(SeqIO.parse(a3m_file,"fasta")) # dict of the modeled fasta, a3m
full_fasta=str(record_dict['1'].seq)

print(full_fasta)
ppdb = PandasPdb()
ppdb.read_pdb(pdb_to_read)

merged_tab=pd.read_table(fasta_tab,sep='\t',index_col=False)

print("i_pdb_id "+i_pdb_id)
print(merged_tab.query('pdb_id == @i_pdb_id')['peptide_fasta'])
pep_fasta_len = len(merged_tab.query('pdb_id == @i_pdb_id')['peptide_fasta'])
prot_fasta_len = len(merged_tab.query('pdb_id == @i_pdb_id')['protein_fasta'])

print("pep_fasta_len "+str(pep_fasta_len))
print("prot_fasta_len "+str(prot_fasta_len))

# figure out which part of the sequence is the peptide and which is the receptor
if pep_fasta_len >1:
    for k in range(pep_fasta_len):
        tmp_seq = str(merged_tab.query('pdb_id == @i_pdb_id')['peptide_fasta'].iloc[k]).upper()
        print("tmp_seq "+tmp_seq)
        if str(tmp_seq) in full_fasta:
            pep_seq = merged_tab.query('pdb_id == @i_pdb_id')['peptide_fasta'].iloc[k]
            break
elif pep_fasta_len ==1:
    pep_seq = merged_tab.query('pdb_id == @i_pdb_id')['peptide_fasta'].iloc[0]
elif pep_fasta_len ==0:
    print("MISSING FASTA")
    seqres_dict = SeqIO.to_dict(SeqIO.parse('/vol/ek/share/peptide_docking_with_afold2_and_rosettAfold/data/pdb_seqres_05022021.txt',"fasta"))
    first_chain_seq = str(seqres_dict[i_pdb_id+"_"+chains[0]].seq)
    second_chain_seq = str(seqres_dict[i_pdb_id+"_"+chains[1]].seq)
    pep_seq = min([first_chain_seq, second_chain_seq], key=len)
    seqres_prot_seq = max([first_chain_seq, second_chain_seq], key=len)
else:
    print("no peptide fasta")

if prot_fasta_len ==0:
    prot_seq=seqres_prot_seq
else:
    prot_seq=str(merged_tab.query('pdb_id == @i_pdb_id')['protein_fasta'].iloc[0]).upper()


if len(pep_seq) > len(prot_seq):
    pep_seq, prot_seq = prot_seq, pep_seq #replace the fastas if for some reason peptide is first
linker_lst=[]

if 'X' in pep_seq:
    print('Xs in pep sequence')
    pep_seq = pep_seq.strip('X')

if 'X' in prot_seq:
    print('Xs in prot sequence')
    prot_seq = prot_seq.strip('X')
    prot_seq = prot_seq.replace('X','')

if int(polyA)!=0:
    pep_seq = len(pep_seq)*'A'

# find the exact positions in the sequence where receptor, linker and peptide begin and end
print("these are full_fasta, prot seq and pep seq")
print(full_fasta)
print(prot_seq)
print(pep_seq)
print("####")
for match in re.finditer(prot_seq,full_fasta): # first and last positions of the receptor
    print('prot match start :', match.start())
    linker_lst.append(match.start())
    print('prot match end :', match.end())
    linker_lst.append(match.end())
    print('linker_lst after prot: ',linker_lst)

for match in re.finditer(pep_seq,full_fasta): # first and last positions of the peptide
    print('pep match start :', match.start())
    print('pep match end :', match.end())
    print(type(match.start()))
    if match.start() >= linker_lst[0] and match.start() <= linker_lst[1]:
        print('match embedded in receptor seq; this is not the peptide')
    else:
        linker_lst.append(match.start())
        linker_lst.append(match.end())

    print('linker_lst after prot and pep: ',linker_lst)

# make a list of all these positions
print('linker_lst includes: ', linker_lst)
df=ppdb.df['ATOM']
linker_lst.append(df[(df['residue_number']<linker_lst[2]) & (df['residue_number']>linker_lst[1])].index[0]) # first linker position, between rec and pep
linker_lst.append(df[(df['residue_number']<linker_lst[2]) & (df['residue_number']>linker_lst[1])].index[-1]) # second linker position, betweem rec and pep

# linker_lst should look like - 
	# receptor_start, receptor_end, pep_start, pep_end, linker_start (index), linker_end (index)

# change the receptor to chain A
receptor_start = df[(df['residue_number']>=linker_lst[0]) & (df['residue_number']<=linker_lst[1])].index[0]
receptor_end = df[(df['residue_number']>=linker_lst[0]) & (df['residue_number']<=linker_lst[1])].index[-1]

df.loc[receptor_start:receptor_start+1,'chain_id']='A'

# change the peptide to chain B
pep_start = df[(df['residue_number']>=linker_lst[2]) & (df['residue_number']<=linker_lst[3])].index[0]
pep_end = df[(df['residue_number']>=linker_lst[2]) & (df['residue_number']<=linker_lst[3])].index[-1]

df.loc[pep_start:(pep_end+1),'chain_id']='B'

# remove the linker
df.drop(df[(df['residue_number']<=linker_lst[2]) & (df['residue_number']>linker_lst[1])].index,axis=0,inplace=True)

# save modified pdb out
ppdb.df['ATOM']=df
ppdb.to_pdb('linker_removed_AB_chains_'+str(pdb_to_read))

