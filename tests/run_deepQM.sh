#! /usr/bin/env bash


deepQM_DIR="/cta/users/otayfuroglu/workspace/deepQM"
PYTHON_DIR="$HOME/miniconda3/bin"


#- set calculation
# Mode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol, opt_single_mol, opt_multi_mol, opt_freq_single_mol
calcMode=opt_grouped_multi_mol

# define model/s: Available model_list="ani1x ani1ccx ani2x aimnetgas aimnetsmd"
model_list="ani2x"


# struct_dir  # pdb files directory
struct_dir=2E


# fbase  for SP calculation
fbase=2E


# namebase 
namebase=2ze1_complex


# set sequence of start (default=0) pose sarting number
# set sequence of end (default=1E6) pose ending number
seq_start=0
seq_end=1


# groups index file
index_file_path="./2E/2ze1_complex.ndx"

# set group1 M23\ protein keyword in pdb
# set group2 (default=SOL) lig or sol keyword in pdb
group1="1"
group2="13"


# set thrshold fmax for optimization (default=0.01)
thr_fmax=1000.01

$PYTHON_DIR/python $deepQM_DIR/deepQM.py $calcMode $model_list $struct_dir $fbase $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax




