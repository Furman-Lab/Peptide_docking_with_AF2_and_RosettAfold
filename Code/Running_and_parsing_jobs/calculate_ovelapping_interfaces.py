#!/usr/bin/env python
# coding: utf-8

# This script calculates the percentage of overlapping interface residues between a native structure and models.

import os
import sys
from pymol import cmd, stored
import numpy as np
import pandas as pd

def get_elements_of_selection(sele):
	myspace = {'res_numbers': []}
	cmd.iterate(sele, 'res_numbers.append(resi)', space=myspace) # iterate through model receptor interface
	output=np.unique(myspace['res_numbers'])

	return(output)

dir = sys.argv[1]
native_file = sys.argv[2]
outdir = sys.argv[3]
outname=sys.argv[4]
all_models = sys.argv[5:]
pdb_id=os.path.basename(os.path.dirname(dir))
complex = pdb_id.split('_')

cmd.set('pdb_no_end_record', 1)
cmd.load(native_file, 'native')
cmd.select("native_rec", 'chain A and backbone')
cmd.select("native_pep", 'chain B and backbone')
cmd.color('white', 'all')

pep_residues = get_elements_of_selection('native_pep')
pep_residues = list(map(int, pep_residues))
first_pos_pep = min(pep_residues)
cmd.alter('native_pep', f'resi=int(resi)-{first_pos_pep}+1')
cmd.select("native_pep", 'chain B and backbone')

# get native interface on peptide and receptor
cmd.select('interface_native_rec', 'native_rec within 8.0 of native_pep and backbone')
native_interface_residues = get_elements_of_selection('interface_native_rec')
native_interface_residues_count = native_interface_residues.size
native_pep_interface_residues = get_elements_of_selection('native_pep within 8 of interface_native_rec and backbone')
native_pep_interface_residues_count = native_pep_interface_residues.size

print('native interface residues' + str(native_interface_residues))
print('native peptide interface residues' + str(native_pep_interface_residues))

metrics = []
for model_file in all_models:
	model_name = os.path.splitext(os.path.basename(model_file))[0]
	cmd.load(model_file, model_name)
	cmd.color('white', model_name)

	cmd.select("afold_rec", f'{model_name} and chain A and backbone')
	cmd.select("afold_pep", f'{model_name} and chain B and backbone')

	cmd.remove('hydrogens')
	cmd.select('interface_afold_rec', 'afold_rec within 8 of afold_pep') # select receptor interface


	# renumber receptor residues
	receptor_residues = get_elements_of_selection('afold_rec')
	receptor_residues = list(map(int, receptor_residues))
	first_pos_receptor = min(receptor_residues)
	cmd.alter('afold_rec', f'resi=int(resi)-{first_pos_receptor}+1')
	cmd.select("afold_rec", f'{model_name} and chain A and backbone')

	# renumber peptide residues
	peptide_residues = get_elements_of_selection('afold_pep')
	peptide_residues = list(map(int, peptide_residues))
	first_pos_peptide = min(peptide_residues)
	cmd.alter('afold_pep', f'resi=int(resi)-{first_pos_peptide}+1')
	cmd.select("afold_pep", f'{model_name} and chain B and backbone')

	afold_interface_residues = get_elements_of_selection('interface_afold_rec')
	afold_interface_residues = list(map(int, afold_interface_residues))
	afold_interface_residues[:] = [ str(number - first_pos_receptor + 1) for number in afold_interface_residues ]
	afold_interface_residues = np.array(afold_interface_residues)
	print('afold interface residues' + str(afold_interface_residues))

	common_residues = np.intersect1d(native_interface_residues, afold_interface_residues) # intersect native and model interface
	common_residues_count = common_residues.size # set of overlapping residues
	common_residues_percent = round(common_residues_count/native_interface_residues_count, 2) # percent of overlapping residues
	print('common_residues', str(common_residues))

	if(common_residues_count > 0):
		print(f'afold_rec and resi ' + '+'.join(common_residues))
		cmd.select('overlapping_interface', f'afold_rec and resi ' + '+'.join(common_residues))
		cmd.select('pep_interact_overlapping_interface', f'afold_pep within 8 of overlapping_interface')

		cmd.color('yellow', 'interface_afold_rec')
		pep_overlapping_interface_residues = get_elements_of_selection('pep_interact_overlapping_interface')
		pep_common_residues = np.intersect1d(pep_overlapping_interface_residues, native_pep_interface_residues)

		print('pep_overlapping_interface_residues', str(pep_overlapping_interface_residues))
		print('pep_common_residues', str(pep_common_residues))
		pep_overlapping_interface_residues_count = pep_common_residues.size # needed to change
		pep_overlapping_interface_residues_percent = round(pep_overlapping_interface_residues_count/native_pep_interface_residues.size, 2)
	else:
		pep_overlapping_interface_residues_count = 0
		pep_overlapping_interface_residues_percent = 0

	metrics.append([pdb_id, model_name, common_residues_count, common_residues_percent,
	pep_overlapping_interface_residues_count, pep_overlapping_interface_residues_percent,
	native_interface_residues_count, native_pep_interface_residues_count])

colnames = ['complex', 'model_name', 'common_residues_count', 'common_residues_percent',
'pep_overlapping_interface_residues_count', 'pep_overlapping_interface_residues_percent',
'native_interface_residues_count', 'native_pep_interface_residues_count']
metrics_df = pd.DataFrame(metrics, columns = colnames)

outfile = os.path.join(outdir, f'{outname}_overlapping_interface.csv')
print(outfile)
metrics_df.to_csv(outfile)
