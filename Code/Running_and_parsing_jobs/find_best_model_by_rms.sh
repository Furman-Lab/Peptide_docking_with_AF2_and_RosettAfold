#!/usr/bin/bash

# script to find the best model based on the rosetta score file and a selected rmsd type

pdb_id=$1

rms_col="number_of_the_relevant_rms_column_in_the_score_file"

# get the best model by the chosen rms column
best_model=$(cat *.sc | grep "my_model" | sort -nk $rms_col | awk '{print $2}' | sed 's/_0001/.pdb/g' | head -1)
cp $best_model best_model.pdb

# grep out only the peptide part (chain B) to compute dssp on - modify to include the relevant model names. e.g. "my_model" ###
grep ' B ' $(cat *.sc | grep "my_model" | sort -nk $rms_col | awk '{print $2}' | sed 's/_0001/.pdb/g' | head -1) > only_peptide_best_model.pdb

