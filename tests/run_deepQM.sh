#! /usr/bin/env bash

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

deepQM_DIR="/cta/users/otayfuroglu/workspace/deepQM"
ANI_DIR="/cta/users/akocak/ASE_ANI/"
PYTHON_DIR="$HOME/miniconda3/bin"

model_list="ani1x ani1ccx ani2x"
$PYTHON_DIR/python $deepQM_DIR/aniQM.py\
	-aniDIR $ANI_DIR\
	-calcMode sp_grouped_multi_mol\
	-model_list $model_list\
	-struct_dir PRO\
	-fbase PRO\
    	-namebase test_\
    	-seq_start 0\
    	-seq_end 8\
    	-group1 ATOM\
    	-group2 HETATM\
    	#-thr_fmax 0.01\
