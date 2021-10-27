#!/usr/bin/bash

### 
# This script gets a directory of models for a specific complex, proccesses the outputs, compares to native and prepares for analyses #
# e.g. 1ssh_AB dir for 1ssh pdb complex, AB = receptor chain, peptide chain; will include the AF2 output for this complex (models, msa files) #



# path to table with modeled pdbs and their fastas (of receptor and peptides)
merged_fasta_tab="path_to_pdb_and_sequence_table"

# define the path to all code
SCRIPTS="path_to_script_directory"
PATH_TO_ROSETTA="path_to_rosetta"
# define the directory being parsed
wd=$1

cd $wd
echo $wd

result=${PWD##*/}

# get the pdb id 
pdb_id=`echo $result | awk -F '_' '{print $1}'`

# get the pdb chains
chains=`echo $result | awk -F '_' '{print $2}'`

# define - is it a polyA run?
polyA=`pwd | grep polyA | wc -l`

msa_file_output_from_af2="name_of_a3m_file_output_from_af2_run.a3m"

echo $pdb_id

# parse models, remove linkers and rename peptide chain
for model in `ls *model*pdb`
do
	if [[ ! $model == *"linker_removed"* ]] # only if was not run as separate chains
	then
		python3 $SCRIPTS/remove_linkers_and_rename_chains.py $model $merged_fasta_tab $msa_file_output_from_af2 $pdb_id $chains $polyA
	fi
done

# get and clean native
native=$pdb_id\_native.pdb

if test -f "$native" # if native already exists, skip this
then
	echo "native exists"
else
	wget -O $pdb_id\_native.pdb 'https://files.rcsb.org/download/'$pdb_id'.pdb'
fi

$PATH_TO_ROSETTA/main/tools/protein_tools/scripts/clean_pdb.py $pdb_id\_native.pdb $chains

# swap the chains if necessary
python3 $SCRIPTS/swap_chains_new.py $pdb_id\_native_$chains\.pdb

# process polyA runs if relevant; for condition to be TRUE and enter the loop, dir path needs to include "polyA" somewhere 
if [[ `pwd | grep polyA | wc -l` -gt 0 ]]
then
	echo "running polyA fixbb on native"
	native_id=`ls *reordered_renamed_native*`
	mv *reordered_renamed_native* original_seq_parsed_native.pdb
	polyA_resfile="path_to_polyA_resfile" # should include a NATRO command and change any AA in the peptide to ALA
	$PATH_TO_ROSETTA/main/source/bin/fixbb.linuxgccrelease -s original_seq_parsed_native.pdb -resfile $polyA_resfile -nstruct 1 -scorefile polyA_score.sc
	mv original_seq_parsed_native_0001.pdb $native_id
fi


native=$(ls *reorder*)
python3 $SCRIPTS/pymol_script.py $PWD ${native} # superimpose to native
python3 $SCRIPTS/remove_unresolved.py $PWD ${native} # remove residues from the models that are unresolved in the native
bash $SCRIPTS/run_flexpepdock.sh $PWD # score and calculdate RMSD with FlexPepDock


# find the best model by rmsd 
bash $SCRIPTS/find_best_model_by_rms.sh $pdb_id

# compute ss2 for models using dssp
python3 $SCRIPTS/compute_dssp_for_models.py

#run alascan on provided model and on native 
bash $SCRIPTS/run_and_parse_alascan_of_best_model_and_native.sh

truncated_pdbs=$(ls *truncated.pdb)

# calculate overlapping interface residues between model and native
	if [[ `echo $truncated_pdbs | wc -l` -gt 0 ]]
        then
                truncated_pdbs=$(echo $truncated_pdbs | tr '\n' ' ')
                python3 $SCRIPTS/calculate_ovelapping_interfaces.py $PWD $PWD/*reorder* $PWD $pdb_id ${truncated_pdbs}
                echo $cmd
                #sbatch --wrap="$cmd"
        fi

