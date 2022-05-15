# Welcome to ANI_LIE calculations!
ANI_LIE is a new tool based on ANI-ML potentials applied on MD trajectories aiming to perform end-state free energy calculations. It has been tested with GROMACS 2018+ versions along with Atomic Simulation Environment.

Please see the documentation below

# Cite us
An accurate binding free energy method from end state MD simulations (just submitted)

The performance of ANI-ML potentials for ligand n(H2O) interaction energies and estimation of hydration free energies from end-point MD simulations (under revision)

# Requirements:
Conda with Python 3
ASE
PyTorch (>=0.4.1)
Torchani
Note: We recommend using a GPU for training the neural networks.

# Installation

conda create --name ANILIE python=3.8

conda activate ANILIE

conda install -c conda-forge numpy pandas tqdm ase pytorch torchani dftd3-python

conda install -c psi4 dftd3

git clone https://github.com/otayfuroglu/deepQM.git

deepQ/tests/run_deepQM.sh 
  
# How To Use

## Tutorial-1: Ligand solvation free energy from ligand+water (LS) simulations

## Tutorial-2: Ligand binding free energy from Protein+ligand+water (PLS) simulations

Although it has been shown that the solvent effects bring little improvements to the accuracy of the calculations, we will provide tutorials for both methods.

### a) Ignoring solvent effects

After MD simulations for protein-ligand aqueous complex (PLS), MD frames can be extracted by using Gromacs trjconv with -sep command in to seperate files. An example trajectory and index file can be found [here](). It is better to create a directory for these frames. 


#### i) ANI only

#### ii) ANI-D3

### b) Including solvent effects

#### i) ANI only

#### ii) ANI-D3
  
  

  
  

  



 
 

