#! /usr/bin/env bash

deepQM_DIR="/cta/users/otayfuroglu/workspace/deepQM"
PYTHON_DIR="/cta/users/eakkus/miniconda3/envs/automd/bin"

#- set calculation
# Mode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol,
# opt_grouped_multi_mol, opt_single_mol, opt_multi_mol
calcMode=sp_grouped_multi_mol

# set number of processors
n_procs=30

# define model/s: Available model_list="ani1x ani1ccx ani2x aimnetgas aimnetsmd"
model_list="ani2x dftd3"

# struct_dir  # pdb files directory
### struct_dir=$(pwd)/pdb
struct_dir=$(pwd)/pdb_pro_lig

# namebase 
namebase=trjmol

# set sequence of start (default=0) pose sarting number
# set sequence of end (default=1E6) pose ending number
seq_start=0
seq_end=-1
###-1 process all pdb files in pdb directory :)

# groups index file
index_file_path="index.ndx"

# set group1 M23\ protein keyword in pdb
# set group2 (default=SOL) lig or sol keyword in pdb
group1="1"
group2="13"

# set thrshold fmax for optimization (default=0.01)
thr_fmax=1000.01
### thr_fmax=0.005

#maximum iteration for optimization
maxiter=5000

$PYTHON_DIR/python $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter
