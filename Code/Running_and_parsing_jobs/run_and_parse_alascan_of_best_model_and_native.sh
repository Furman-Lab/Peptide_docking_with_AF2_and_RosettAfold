#!/usr/bin/bash

# script to be run from within a directory that is being processed. creates another directory witihin the working dir, copies the relevant model and native,
# runs rosetta alascan on both

mkdir alascan_best_model_and_native
cd alascan_best_model_and_native
alascan="path_to_rosetta_alascan_bin"

# A and B stand for renamed chains of the receptor and peptide respectively. change if other chains are used
$alascan ../best_model A B
$alascan ../*reordered_renamed_native.pdb A B
cd ../


