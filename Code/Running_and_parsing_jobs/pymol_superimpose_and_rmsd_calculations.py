#!/bin/python3

# This script superimposes the models by their receptor and calculates the receptor and
# peptide RMSD both by superimposing and aligning by sequence. It also creates a PyMol
# session with the interface residues colored. Input id a director with a native and
# models that start with 'linker_removed'.

from pymol import cmd, stored
import interfaceResidues
import pandas as pd
import os
import re
import glob
import Focus_alignment
import sys

# load native and select receptor and peptide
input_dir=sys.argv[1]
print(input_dir)
complex=os.path.basename(input_dir.rstrip('/'))
splitted_complex = complex.split('_', 1)
pdb_id = splitted_complex[0]
rec_chain = 'A'
pep_chain = 'B'
cmd.load(sys.argv[2], 'native')

# remove not receptor and peptide chains
cmd.select('to_remove', f'native and not chain {rec_chain} and not chain {pep_chain}')
cmd.remove('to_remove')
cmd.remove('resn HOH')
cmd.remove('hydrogens')

# select peptide and receptor chain
cmd.select("native_rec", f'native and chain {rec_chain}')
cmd.select("native_pep", f'native and chain {pep_chain}')

# select interface by first selecting receptor residues within 4A of peptide, then selecting peptide residues within 4A of receptor interface
cmd.select('native_rec_interface', 'byres native_rec within 4 of native_pep')
cmd.select('interface_native', 'native_rec_interface + byres native_pep within 4 of native_rec_interface')

# color receptor interface of native in yellow
cmd.select('interface_rec', 'interface_native and chain {rec_chain}')
cmd.color('orange', 'interface_rec')

# color peptide interface of model in cyan
cmd.color('red', 'native_pep')
cmd.show(representation='sticks', selection='native_pep')


# setting up common variables
metrics = []
rank_re = re.compile("rank_[0-9]")
model_re = re.compile("model_[0-9]")


# load models 1by1 select receptor and peptide
for model_file in glob.glob(os.path.join(input_dir, 'linker_removed*.pdb')):
	model_name = os.path.basename(model_file).replace('.pdb', '')
	rank = int(str(rank_re.search(model_name).group(0)).replace('rank_', ''))
	model_no = int(str(model_re.search(model_name).group(0)).replace('model_', ''))
	metrics_for_a_model = [model_name, complex, rec_chain, pep_chain, rank, model_no]

	cmd.load(model_file, model_name)
	cmd.select("afold_rec", f'{model_name} and chain A')
	cmd.select("afold_pep", f'{model_name} and chain B')

	# align peptide chains
	super_alignment_pep = cmd.super('afold_pep', 'native_pep')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_pep])

	seq_alignment_pep = cmd.align('afold_pep', 'native_pep')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_pep])

	# align peptide chains backbone
	super_alignment_pep = cmd.super('afold_pep and backbone', 'native_pep and backbone')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_pep])

	seq_alignment_pep = cmd.align('afold_pep and backbone', 'native_pep and backbone')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_pep])


	# super receptor chains
	super_alignment_rec = cmd.super('afold_rec', 'native_rec')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_rec])

	# save the superimposed structure
	cmd.select('model_to_save', model_name)
	super_filename=f'{model_name}_superimposed.pdb'
	print(super_filename)
	save_to_file = os.path.join(input_dir, super_filename)
	cmd.save(save_to_file, model_name, format='pdb')

	# super receptor chain backbones
	super_alignment_rec = cmd.super('afold_rec and backbone', 'native_rec and backbone')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_alignment_rec])

	# calculate rmsd-s
	seq_alignment_rec = cmd.align('afold_rec', 'native_rec')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_rec])

	# calculate rmsd by backbone
	seq_alignment_rec = cmd.align('afold_rec and backbone', 'native_rec and backbone')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in seq_alignment_rec])

	super_complex = cmd.super(model_name, 'native')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_complex])

	super_complex = cmd.super(f'{model_name} and backbone', 'native and backbone')
	metrics_for_a_model += tuple([float("{0:.2f}".format(n)) for n in super_complex])

	cmd.color('brown', model_name)
	# color receptor interface of model in yellow
	cmd.color('yellow', 'interface_rec_afold')

	# color peptide of model in cyan
	cmd.color('cyan', 'afold_pep')
	cmd.show(representation='sticks', selection='afold_pep')

	metrics.append(metrics_for_a_model)

cmd.set_name('native', complex)
cmd.save(f'{input_dir}/{complex}.pse', format='pse')

# create column names
colnames = ['model_name', 'pdb_id', 'rec_chain', 'pep_chain', 'rank', 'model_no']
colnames_for_aln = ['rms_after_ref', 'no_aln_atoms_after_ref', 'ref_cycles', 'rms_before_ref', 'no_aln_atoms_before_ref', 'raw_score', 'no_aln_residues']

for type in ['_super_pep', '_align_seq_pep', '_super_pep_bb', '_align_seq_pep_bb', '_super_rec', '_super_rec_bb', '_align_seq_rec', '_align_seq_rec_bb', '_complex', '_complex_bb']:
	new_colnames = [s + type for s in colnames_for_aln]
	colnames = colnames + new_colnames


# saving calculated metrics
out_csv_dir = os.path.dirname(input_dir.rstrip('/'))
metrics_df = pd.DataFrame(metrics, columns = colnames)
metrics_df.to_csv(os.path.join(out_csv_dir, complex + '.csv'))
