#! /usr/bin/env bash
# You can run the script on CUDA or CPU
export CUDA_VISIBLE_DEVICES=0  # one GPU usage
# export CUDA_VISIBLE_DEVICES=""  # run on CPU

# You need to set deepQM and python path for your case
deepQM_DIR="$HOME/deepQM"
PYTHON_DIR="$HOME/miniconda3/envs/ANILIE/bin"

#- set calculation type
# Mode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol,
# opt_grouped_multi_mol, opt_single_mol, opt_multi_mol
calcMode=sp_grouped_multi_mol

# set number of parallel tasks. It can be GPU or CPU. If using GPU, the GPU-RAM may be limitatition so reduce n_procs for such cases.
n_procs=2

# define model/s: Available model_list=ani1x, ani1ccx, ani2x, g16, dftd3. One or more models can be set at once.
model_list="ani2x dftd3"

# The location of pdb files. Explicit directory or subdirectories can be defined.
struct_dir=$(pwd)/pdb_pro_lig_sol

# This tells what the base name of the pdb files are.
# When the trajectory is extracted in to multiple files, the file name excluding the frame number is given here. 
# (e.g. trjmol1.pdb, trjmol2.pdb,... basename=trjmol)
namebase=pls

# frame starting number for the basename (default=0).
seq_start=0

# frame ending number (default=1000). if -1 is given, it calculates all the pdb files in the struct_dir without caring namebase.
seq_end=1001

# Group index file in Gromacs format. 
# Any two groups can be chosen from index file. 
# Required when "group" types of calculators are set. 
# Groups are defined in the following lines
# IMPORTANT: Full path must be given.
index_file_path=$(pwd)/index_4A.ndx
group1="13"
group2="15"

# sets threshold fmax for optimization (default=0.01). Only used in optimization types of calculators
thr_fmax=0.01

# Maximum iteration for optimization. Only used in optimization types of calculators
maxiter=5000

#$PYTHON_DIR/python $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter

# For ANI_LIE: alpha=0.000, beta=0.272 and gamma=-2.164. For ANID3_LIE: alpha=-0.057, beta=0.208 and gamma=-1.230.
$PYTHON_DIR/python $deepQM_DIR/scripts/bindEnAniD3.py -in "$struct_dir"/"$namebase"_SP_energies_"$group1"_"$group2".csv -ani ani2x -a -0.057 -b 0.208 -g -1.230
