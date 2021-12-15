#! /usr/bin/env bash

# set ani environment variables
export PATH=/usr/local/cuda-9.2/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-9.2/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

###### NC_ROOT EXPORTS ########
export NC_ROOT="/home/modellab/ASE_ANI"
export LD_LIBRARY_PATH="$NC_ROOT/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="$NC_ROOT/lib:$PYTHONPATH"

# keyword list:
#-calcMode: Options -> sp_single_mol, sp_multi_mol, sp_grouped_multi_mol, opt_single_mol, opt_multi_mol, opt_freq_single_mol
#-struct_dir  # pdb files directory
#-fbase  for SP calculation
#-namebase 
#-seq_start (default=0) pose sarting number
#-seq_end (default=1E6) pose ending number
#-group1 M23\ protein keyword in pdb
#-group2 (default=SOL) lig or sol keyword in pdb
#-thr_fmax (default=0.01)

ANI_DIR="/home/modellab/ASE_ANI"
PYTHON_DIR="$HOME/miniconda3/envs/automd/bin"
$PYTHON_DIR/python aniQM.py\
	-aniDIR $ANI_DIR\
	-calcMode sp_grouped_multi_mol\
	-struct_dir /home/modellab/Desktop/PROJE_serdar/file\
	-fbase .\
    	-namebase deneme_\
    	-seq_start 1\
    	-seq_end 2\
    	-group1 ATOM\
    	-group2 HETATM\
    	#-thr_fmax 0.01\
