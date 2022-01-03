#! /usr/bin/env bash


deepQM_DIR="/cta/users/otayfuroglu/workspace/deepQM"
PYTHON_DIR="$HOME/miniconda3/bin"


#- set calculation
# Mode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol,
# opt_grouped_multi_mol, opt_single_mol, opt_multi_mol
calcMode=sp_grouped_multi_mol

# set number of processors
n_procs=1

# define model/s: Available model_list="ani1x ani1ccx ani2x aimnetgas aimnetsmd"
model_list="ani1x ani1ccx ani2x"


# struct_dir  # pdb files directory
struct_dir=2E


# namebase 
namebase=2ze1_complex


# set sequence of start (default=0) pose sarting number
# set sequence of end (default=1E6) pose ending number
seq_start=1
seq_end=2


# groups index file
index_file_path="$struct_dir/index.ndx"
# index_file_path="None"

# set group1 M23\ protein keyword in pdb
# set group2 (default=SOL) lig or sol keyword in pdb
group1="1"
group2="13"


# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.7

#maximum iteration for optimization
maxiter=500

$PYTHON_DIR/python -W ignore $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter




