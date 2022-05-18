#! /usr/bin/env bash


deepQM_DIR="/cta/users/otayfuroglu/workspace/deepQM"
PYTHON_DIR="$HOME/miniconda3/bin"


#- set calculation
# Mode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol,
# opt_grouped_multi_mol, opt_single_mol, opt_multi_mol
calcMode=sp_grouped_multi_mol

# set number of processors. It can be GPU or CPU. If using GPU, the GPU-RAM may be limitatition so reduce n_procs for such cases.
n_procs=2

# define model/s: Available model_list="ani1x ani1ccx ani2x aimnetgas aimnetsmd"
model_list="ani2x dftd3"


# struct_dir  # includes pdb files dir.
struct_dir=/cta/users/otayfuroglu/workspace/deepQM/tests/pdb_pro_lig


# namebase 
namebase=trjmol


# set sequence of start (default=0) pose sarting number
# set sequence of end (default=1001) pose ending number
seq_start=0
seq_end=1001


# groups index file
index_file_path=$struct_dir/index.ndx
# index_file_path="None"

# set group1  (default=Protein) keyword in index
# set group2 (default=SOL) lig or sol keyword in index
group1="1"
group2="13"


# set thrshold fmax for optimization (default=0.01)
thr_fmax=0.05

#maximum iteration for optimization
maxiter=500

$PYTHON_DIR/python -W ignore $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter
$PYTHON_DIR/python $deepQM_DIR/scripts/bindEnAniD3.py -in "$struct_dir"/"$namebase"_SP_energies_"$group1"_"$group2".csv -ani ani2x -a 0.0 -b 0.127 -g -5.111



