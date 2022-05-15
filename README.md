# Welcome to ANI_LIE calculations!
ANI_LIE is a new tool based on ANI-ML potentials applied on MD trajectories aiming to perform end-state free energy calculations. It has been tested with GROMACS 2018+ versions along with Atomic Simulation Environment.

Please see the documentation below

# Cite us
An accurate binding free energy method from end state MD simulations (just submitted)

The performance of ANI-ML potentials for ligand n(H2O) interaction energies and estimation of hydration free energies from end-point MD simulations (under revision)

# Requirements:
- Conda with Python 3
- ASE
- PyTorch (>=0.4.1)
- Torchani
Note: We recommend using a GPU for training the neural networks.

# Installation
```
conda create --name ANILIE python=3.8

conda activate ANILIE

conda install -c conda-forge numpy pandas tqdm ase pytorch torchani dftd3-python

conda install -c psi4 dftd3

git clone https://github.com/otayfuroglu/deepQM.git

cd deepQ/tests/run_deepQM.sh 
```
# How To Use

deepQM.py is the main pyhton script that does the calculations. It can be called from run_deepQM.sh file in which a set of parameters are predefined and can ben customized according to the calculation type.

## Parameters
### calcMode
There are several types of calculation modes: sp_single_mol, sp_multi_mol, sp_grouped_multi_mol, opt_grouped_multi_mol, opt_single_mol, opt_multi_mol

sp_single_mol:single point energy of a single pdb file

sp_multi_mol: single point energy of a multiple pdb files. It calls sp_single_mol for all compounds in a directory

sp_grouped_multi_mol: single point energy of a multiple compounds which runs the sp_single_mol three times and finds differences between two groups (group A and group B) and it performs the same calculation for multiple pdb files

opt_single_mol: it performs geometry optimization for a single pdb file

opt_multi_mol: it performs geometry optimization for multiple pdb files

### model_list
This sets the calculator type currently ani1x, ani1ccx, ani2x, dftd3 and g16 calculators are provided. More than one type of calculators can be set. The resulting csv file will print all the modellist.
### struct_dir
This tells where the pdb files are. Explicit directory or subdirectories can be defined.
### namebase
This tells what the base name of the pdb files. When the trajectory is extracted in to multiple files the file name excluding the frame number is given here.
This name is only required when sequential file names are given with seq_start and seq_end commands.
### seq_start
sets sequence of start (default=0) frame starting number. 
### seq_end
sets frame ending number (default=1000). if -1 is given, it calculates all the pdb files in the struct_dir without caring namebase.
### index_file_path
Group index file in Gromacs format. Any two groups can be chosen from index file. Required when "group" types of calculators. Groups are defined in the following lines
### group1
Group A in AB-->A+B calculation Default is 1 for protein in MD simulations
### group2
Group B in AB-->A+B calculation Default is 13 for ligand in MD simulations
### thr_fmax
set thrshold fmax for optimization (default=0.01). Only used in optimization types of calculators
### maxiter
Maximum iteration for optimization. Only used in optimization types of calculators

### run script by:
```
$PYTHON_DIR/python $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter
```
## Tutorial-1: Ligand solvation free energy from ligand+water (LS) simulations


## Tutorial-2: Ligand binding free energy from Protein+ligand+water (PLS) simulations

Although it has been shown that the solvent effects bring little improvements to the accuracy of the calculations, we will provide tutorials for both methods.

### a) Ignoring solvent effects




#### i) ANI only

After MD simulations for protein-ligand aqueous complex (PLS), MD frames can be extracted by using Gromacs trjconv with -sep command in to seperate files. An example MD simulation with trajectory and index file can be found [here](). It is better to create a directory for these frames. We have to make sure periodic boundary condition is removed from trajectory

```
mkdir pdb_pro_lig
echo 24 24 | gmx trjconv -f md_0.xtc -s md_0.tpr -n index.ndx -o pdb_pro_lig/trjmol.pdb -pbc nojump -ur compact -center -sep
```
Here index group 24 is the Protein_Ligand complex while 1 is the protein and 13 is the ligand. You can use -dt option to reduce number of frames to calculate. Let's say we have 1000 frames (trjmol0.pdb, trjmol1.pdb, ... trjmol1000.pdb) extracted into a directory name of ./pdb_pro_lig/
```
bash run_deepQM.sh [here](tests/run_deepQM.sh)
```
After running this command, it will create a csv file in the same directory of ./pdb_pro_lig/. This file includes all frames for each group and energy differences written in eVs without fitting coefficients.
```
$PYTHON_DIR/python $deepQM_DIR/scripts/bindEnAniD3.py -in "$struct_dir"/"$namebase"_SP_energies_"$group1"_"$group2".csv -ani ani2x -a 0.0 -b 0.127 -g -5.111
```
This script collects the data from the csv file previously produced and converts to free energies with coefficients determined from fit to the experimental energies







#### ii) ANI-D3

### b) Including solvent effects

coming soon

#### i) ANI only

coming soon

#### ii) ANI-D3
  
coming soon
  

  
  

  



 
 

