#! /usr/bin/env bash
# You can run the script on CUDA or CPU
export CUDA_VISIBLE_DEVICES=0  # one GPU usage
# export CUDA_VISIBLE_DEVICES=""  # run on CPU


deepQM_DIR="/cta/users/eakkus/deepQM"
PYTHON_DIR="/cta/users/eakkus/miniconda3/envs/automd/bin"


#- set calculation
# Mode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol,
# opt_grouped_multi_mol, opt_single_mol, opt_multi_mol
calcMode=sp_grouped_multi_mol

# set number of processors
n_procs=2

# define model/s: Available model_list="ani1x ani1ccx ani2x aimnetgas aimnetsmd"
model_list="ani2x dftd3"

# struct_dir  # pdb files directory
### struct_dir=$(pwd)/pdb
struct_dir=$(pwd)/pdb_pro_lig_sol

# namebase 
namebase=trjmol

# set sequence of start (default=0) pose sarting number
# set sequence of end (default=1E6) pose ending number
seq_start=0
seq_end=1001
###-1 process all pdb files in pdb directory :)

# groups index file
index_file_path="index_4A.ndx"

group1="13"
group2="15"

# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.01
### thr_fmax=0.005

#maximum iteration for optimization
maxiter=5000

#$PYTHON_DIR/python $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter
$PYTHON_DIR/python $deepQM_DIR/scripts/bindEnAniD3.py -in "$struct_dir"/"$namebase"_SP_energies_"$group1"_"$group2".csv -ani ani2x -a -0.057 -b 0.208 -g -1.230