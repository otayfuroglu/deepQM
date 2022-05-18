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

cp deepQ/tests/run_deepQM.sh <working-dir>
```
# How To Use

deepQM.py is the main pyhton script that does the calculations. It can be called from run_deepQM.sh file in which a set of parameters are predefined and can be customized according to the calculation type.

## Parameters
### calcMode
There are several types of calculation modes: sp_single_mol, sp_multi_mol, sp_grouped_multi_mol, opt_grouped_multi_mol, opt_single_mol, opt_multi_mol. Default is sp_grouped_multi_mol.

sp_single_mol:single point energy of a single pdb file

sp_multi_mol: single point energy of multiple pdb files. It calls sp_single_mol for all pdb files in a directory.

sp_grouped_multi_mol: single point energy of a multiple compounds which runs the sp_single_mol three times and finds differences between two groups (group A and group B) and it performs the same calculation for multiple pdb files

opt_single_mol: geometry optimization for a single pdb file

opt_multi_mol: geometry optimization for multiple pdb files

### model_list
This sets the calculator type. Currently ani1x, ani1ccx, ani2x, dftd3 and g16 calculators are provided. More than one type of calculators can be set. The resulting csv file will print all the model_list calculator results. The default is ani2x
### struct_dir
This tells where the pdb files are. Explicit directory or subdirectories can be defined.
### namebase
This tells what the base name of the pdb files are. When the trajectory is extracted in to multiple files, the file name excluding the frame number is given here. (e.g. trjmol1.pdb, trjmol2.pdb,... basename=trjmol)
This name is only required when sequential file names are given with seq_start and seq_end commands.
### seq_start
sets frame starting number for the basename (default=0). 
### seq_end
sets frame ending number (default=1000). if -1 is given, it calculates all the pdb files in the struct_dir without caring namebase.
### index_file_path
Group index file in Gromacs format. Any two groups can be chosen from index file. Required when "group" types of calculators are set. Groups are defined in the following lines
### group1
Group A in AB-->A+B calculation Default is 1 for protein in MD simulations
### group2
Group B in AB-->A+B calculation Default is 13 for ligand in MD simulations
### thr_fmax
sets thrshold fmax for optimization (default=0.01). Only used in optimization types of calculators
### maxiter
Maximum iteration for optimization. Only used in optimization types of calculators

### run script by:
```
$PYTHON_DIR/python $deepQM_DIR/deepQM.py $calcMode $n_procs $model_list $struct_dir $namebase $seq_start $seq_end $index_file_path $group1 $group2 $thr_fmax $maxiter
```
### analyze results by: 
```
$PYTHON_DIR/python $deepQM_DIR/scripts/bindEnAniD3.py -in "$struct_dir"/"$namebase"_SP_energies_"$group1"_"$group2".csv -ani ani2x -a 0.0 -b 0.127 -g -5.111
```
## Tutorial-1: Ligand solvation free energy from ligand+water (LS) simulations

After MD simulations for ligand aqueous complex (LS), MD frames can be extracted by using Gromacs trjconv with -sep command in to seperate files. An example MD simulation with trajectory and index file can be found [here](tests/LS/). It is better to create a directory for these frames (pdb_lig_sol). We have to make sure periodic boundary condition is removed from trajectory.

```
mkdir pdb_lig_sol
echo 2 0 | gmx trjconv -f md_0.xtc -s md_0.tpr -n index.ndx -o pdb_lig_sol/trjmol.pdb -pbc mol -ur compact -center -sep
```
Here index group 2 is the ligand while 0 is the system. You can use -dt option to reduce number of frames to calculate. We prefer to have at least 100-200 frames. Let's say we have 1000 frames (trjmol0.pdb, trjmol1.pdb, ... trjmol1000.pdb) extracted into a directory name of ./pdb_lig_sol/.
```
bash run_deepQM.sh
```
After running this command, it will create a csv file in the same directory of ./pdb_lig_sol/. This file includes all frames for each group and energy differences written in eVs without fitting coefficients.

The last line of the script collects the data from the csv file previously produced and converts to free energies with coefficients determined from fit to the experimental energies using AVG/ANI_LIE. 

DG=alpha*<diff_dftd3>+beta*<diff_ani2x>+gamma formula is used in the calculations.

For ANI_LIE: alpha=0.000, beta=0.272 and gamma=-2.164. For ANID3_LIE: alpha=-0.057, beta=0.208 and gamma=-1.230. Experimental values are detrmined from fit. The values may change slightly according to differen system.

If the user wants to calculate the other methods discussed elsewhere (EXP, nc-EXP and cu-EXP), we have a separate script that can calculate all different methods [here](scripts/statsdeepAniOutputs.py). 

## Tutorial-2: Ligand binding free energy from Protein+ligand+water (PLS) simulations

Although it has been shown that the solvent effects bring little improvements to the accuracy of the calculations, we will provide tutorials for both methods.

### a) Ignoring solvent effects

After MD simulations for protein-ligand aqueous complex (PLS), MD frames can be extracted by using Gromacs trjconv with -sep command in to seperate files. An example MD simulation with trajectory and index file can be found [here](). It is better to create a directory for these frames. We have to make sure periodic boundary condition is removed from trajectory

```
mkdir pdb_pro_lig
echo 24 24 | gmx trjconv -f md_0.xtc -s md_0.tpr -n index.ndx -o pdb_pro_lig/trjmol.pdb -pbc nojump -ur compact -center -sep
```
Here index group 24 is the Protein_Ligand complex while 1 is the protein and 13 is the ligand. You can use -dt option to reduce number of frames to calculate. We prefer to have at least 100-200 frames. Let's say we have 1000 frames (trjmol0.pdb, trjmol1.pdb, ... trjmol1000.pdb) extracted into a directory name of ./pdb_pro_lig/.
```
bash run_deepQM.sh
```
After running this command, it will create a csv file in the same directory of ./pdb_pro_lig/. This file includes all frames for each group and energy differences written in eVs without fitting coefficients.

The last line of the script collects the data from the csv file previously produced and converts to free energies with coefficients determined from fit to the experimental energies. 

DG=alpha*<diff_dftd3>+beta*<diff_ani2x>+gamma formula is used in the calculations.

For ANI_LIE: alpha=0.000, beta=0.127 and gamma=-5.111. For ANID3_LIE: alpha=-0.0353, beta=0.1487 and gamma=-5.9866. Experimental values are detrmined from fit. The values may change slightly according to differen system.



### b) Including solvent effects

In order to account for solvation terms, we need another MD simulation of free ligand in water (LS simulation) in addition to Protein+ligand+water (PLS) simulation. We will find the binding energy in three steps.

#### step-1) LS simulation for LSinLS

This part is already discussed in [LS simulation](https://github.com/otayfuroglu/deepQM#tutorial-1-ligand-solvation-free-energy-from-ligandwater-ls-simulations). Obtained result, which we call LSinLS from [this script](tests/LS/run_deepQM.sh) will be used in step-3.

#### step-2) PLS simulation for LSinPLS
  
Ideally, after running PLS simulaions, you can calculate L-surr by creating groups of L vs PS from the PLS simulations. Equivalently, you can use L vs P and L vs S. You will need to run the script twice for the P-L and L-S interactions by creating corresponding index groups. The second way requires less memory use in the computations, so we will go with this one. Having all waters in the index groups is still too costly in particular for DFTD3 calculations. Here we will use reduced systems (considering only 4A waters around the PL complex) but if the RAM of your GPU/CPU allows you can use all waters.

To reduce your system size, you can use following Gromacs trjorder command. But be aware that new index groups for these need to be generated for this reduced system. We have a [script](tests/PLS/PLS_4A.sh) that does 
-remove periodicity

-order trajectory

-assign number of solvents for 4A around PL complex

-generate new index group for solvents within 4A

-create trajectory of reduced system



In the first run, you will still run the same script discussed in [LS simulation](https://github.com/otayfuroglu/deepQM#tutorial-1-ligand-solvation-free-energy-from-ligandwater-ls-simulations)

  
  

  
  

  



 
 

