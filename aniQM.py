#
import sys
from ase_interface import ANIENS
from ase_interface import aniensloader

import  ase
import time
from ase import units
from ase.io import read, write, trajectory
from ase.optimize import BFGS, LBFGS
from ase.vibrations import Vibrations

import pandas as pd
import os, re
from operator import itemgetter
import tqdm

import numpy  as np
import argparse

from aniqm_utils import prepare_xyz_grp, prepare_xyz_complex, pdb2xyz

parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-aniDIR", "--aniDIR",
                    type=str, required=True,
                    help="")
parser.add_argument("-calcMode", "--calcMode",
                    type=str, required=True,
                    help="")
parser.add_argument("-struct_dir", "--struct_dir",
                    type=str, required=True,
                    help="")
parser.add_argument("-fbase", "--fbase",
                    type=str, required=True,
                    help="")
parser.add_argument("-namebase", "--namebase",
                    type=str, required=False,
                    help="")
parser.add_argument("-seq_start", "--seq_start",
                    type=int, required=False, default=0,
                    help="")
parser.add_argument("-seq_end", "--seq_end",
                    type=int, required=False, default=1E6,
                    help="")

parser.add_argument("-group1", "--group1",
                    type=str, required=False,
                    help="")

parser.add_argument("-group2", "--group2",
                    type=str, required=False, default="SOL",
                    help="")

parser.add_argument("-thr_fmax", "--thr_fmax",
                    type=float, required=False, default=0.01,
                    help="")


args = parser.parse_args()
# set envs
aniDIR = args.aniDIR
model_anix_path = '%s/ani_models/ani-1x_8x.info' %aniDIR
model_aniccx_path = '%s/ani_models/ani-1ccx_8x.info' %aniDIR
model_ani2x_path = '%s/ani_models/ani-2x_8x.info' %aniDIR

def calcSPWithAni(model_path, traj_path, mol):

    # mol_path = structure_dir + "/" + file_base + ".xyz"

    # check .xyz structure file
    #if not os.path.exists(mol_path):
    #    pdb2xyz(structure_dir, file_base)

    # Read molecule from xyz
    # mol = read(mol_path)
    # set calculator
    mol.set_calculator(ANIENS(aniensloader(model_path, 0)))

    start_time = time.time()

    # Calculate energy
    ef = mol.get_potential_energy()
    run_time = time.time() - start_time
    return float("%.8f"%ef)

def runSPgroupedMultiMol(grp1, grp2, directory, base, seq_start, seq_end):
    list_results = []
    structure_dir = directory
    traj_path =  "%s/traj_files" %structure_dir
    
    if not os.path.exists(traj_path):
        traj_path = os.mkdir(traj_path)

    labels = ["Structure", "SPE_anix_%s_%s"%(grp1, grp2), "SPE_anix_%s"%grp1, "SPE_anix_%s"%grp2,
              "SPE_aniccx_%s_%s"%(grp1, grp2), "SPE_aniccx_%s"%grp1, "SPE_aniccx_%s"%grp2,
              "SPE_ani2x_%s_%s"%(grp1, grp2), "SPE_ani2x_%s"%grp1, "SPE_ani2x_%s"%grp2,
              "diff_anix", "diff_aniccx", "diff_ani2x"]
    df_results = pd.DataFrame(columns=labels)
    df_results.to_csv("%s/%s_SP_energies.csv" %(structure_dir, base))

    for i in range(seq_start, seq_end):
        file_base = "{}{}".format(base, i)

        prepare_xyz_complex(structure_dir, file_base, grp1, grp2)
        mol_path = structure_dir + "/" + file_base + ".xyz"
        mol = read(mol_path)
        atom_types = mol.get_chemical_symbols()

        list_results_anix = [str(file_base)]
        list_results_aniccx = []
        list_results_ani2x = []
        #i += 1
        print("%d. pdb file is processing ..." %i)
        #file_base = file_name.replace(".pdb", "").replace(".xyz", "")
        print("Calcultion starts for %s" %file_base)

        # check atom type for ani2x
        ani2_atom_type = [sym for sym in atom_types if sym in ["Cl", "F", "S"]]

        if len(ani2_atom_type) >= 1:
            results_anix = 0
            results_aniccx = 0
            results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
        else:
            results_anix = calcSPWithAni(model_anix_path, traj_path, mol)
            results_aniccx = calcSPWithAni(model_aniccx_path, traj_path, mol)
            results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
        list_results_anix.append(results_anix)
        list_results_aniccx.append(results_aniccx)
        list_results_ani2x.append(results_ani2x)

        for grp in [grp1, grp2]:
            prepare_xyz_grp(grp, structure_dir, file_base)
            file_base_new = "{}_{}".format(file_base, grp)
            mol_path = structure_dir + "/" + file_base_new + ".xyz"
            mol = read(mol_path)
            if len(ani2_atom_type) >= 1:
                results_anix = 0
                results_aniccx = 0
                results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
            else:
                results_anix = calcSPWithAni(model_anix_path, traj_path, mol)
                results_aniccx = calcSPWithAni(model_aniccx_path, traj_path, mol)
                results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
            list_results_anix.append(results_anix)
            list_results_aniccx.append(results_aniccx)
            list_results_ani2x.append(results_ani2x)
        list_results = list_results_anix + list_results_aniccx + list_results_ani2x\
                + [list_results_anix[1] - list_results_anix[2] - list_results_anix[3]] \
                + [list_results_aniccx[0] - list_results_aniccx[1] - list_results_aniccx[2]] \
                + [list_results_ani2x[0] - list_results_ani2x[1] - list_results_ani2x[2]]
            #print(list_results)
        #  except:
        #      print("Error for {}\nSkipping...\n".format(file_base))
        #      list_results = np.zeros(len(labels))
        #      list_results[0] = str(file_base)
        #print(list_results)
        df_results = pd.DataFrame([list_results], columns=labels)
        df_results.to_csv("%s/%s_SP_energies.csv" %(structure_dir, base), mode="a", header=False)

def runSPmultiMol(directory, base, seq_start, seq_end):
    list_results = []
    structure_dir = directory
    traj_path =  "%s/traj_files" %structure_dir

    if not os.path.exists(traj_path):
        traj_path = os.mkdir(traj_path)

    labels = ["Structure", "SPE_anix", "SPE_aniccx", "SPE_ani2x"]
    df_results = pd.DataFrame(columns=labels)
    df_results.to_csv("%s/%s_SP_energies.csv" %(structure_dir, directory))

    for i in tqdm.trange(seq_start, seq_end):
        file_base = "{}{}".format(base, i)
        mol_path = structure_dir + "/" + file_base + ".xyz"
        mol = read(mol_path)
        atom_types = mol.get_chemical_symbols()
        #i += 1
        #  print("%d. pdb file is processing ..." %i)
        #  print("Calcultion starts for %s" %file_base)

        # check atom type for ani2x
        ani2_atom_type = [sym for sym in atom_types if sym in ["Cl", "F", "S"]]

        if len(ani2_atom_type) >= 1:
            results_anix = 0
            results_aniccx = 0
            results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
            list_results = [str(file_base), results_anix, results_aniccx, results_ani2x]
        else:
            results_anix = calcSPWithAni(model_anix_path, traj_path, mol)
            results_aniccx = calcSPWithAni(model_aniccx_path, traj_path, mol)
            results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
            list_results = [str(file_base), results_anix, results_aniccx, results_ani2x]

        df_results = pd.DataFrame([list_results], columns=labels)
        df_results.to_csv("%s/%s_SP_energies.csv" %(structure_dir, directory), mode="a", header=False)

def runSP(directory, file_base):
    list_results = []
    structure_dir = directory
    traj_path =  "%s/traj_files" %structure_dir

    mol_path = structure_dir + "/" + file_base + ".xyz"
    mol = read(mol_path)
    atom_types = mol.get_chemical_symbols()

    if not os.path.exists(traj_path):
        traj_path = os.mkdir(traj_path)

    # check atom type for ani2x
    ani2_atom_type = [sym for sym in atom_types if sym in ["Cl", "F", "S"]]

    if len(ani2_atom_type) >= 1:
        results_anix = 0
        results_aniccx = 0
        results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
        list_results = [str(file_base), results_anix, results_aniccx, results_ani2x]
    else:
        results_anix = calcSPWithAni(model_anix_path, traj_path, mol)
        results_aniccx = calcSPWithAni(model_aniccx_path, traj_path, mol)
        results_ani2x = calcSPWithAni(model_ani2x_path, traj_path, mol)
        list_results = [str(file_base), results_anix, results_aniccx, results_ani2x]

    labels = ["Mol", "SPE_anix", "SPE_aniccx", "SPE_ani2x"]
    for label, value in zip(labels, list_results):
        print("{} --> {}".format(label, value))

def runOptSingleMol(structure_dir, file_base, thr_fmax):

    traj_path =  "%s/traj_files" %structure_dir

    if not os.path.exists(traj_path):
        traj_path = os.mkdir(traj_path)
    
    prepare_xyz_complex(structure_dir, file_base)

    mol_path = structure_dir + "/" + file_base + ".xyz"

    # check .xyz structure file
    #if not os.path.exists(mol_path):
    #pdb2xyz(structure_dir, file_base)
    
    # Read molecule from xyz
    mol = read(mol_path)

    mol.set_calculator(ANIENS(aniensloader(model_ani2x_path, 0)))

    start_time = time.time()
    #dyn = LBFGS(mol, trajectory="%s/%s.traj" %(traj_path, file_base))
    dyn = LBFGS(mol)
    dyn.run(fmax=thr_fmax)
    run_time = time.time() - start_time
    print('[ANI Optimization - Total time:', run_time, 'seconds]')

    # Calculate energy
    ef = mol.get_potential_energy()
    write("%s/optimzed_%s.xyz" %(structure_dir, file_base), mol)
    print("Final Energy: ", ef)
    print("job is end for %s" %file_base)
    return [ef, run_time]

def runOptFreqSingleMol(structure_dir, file_base, thr_fmax):

    traj_path =  "%s/traj_files" %structure_dir

    if not os.path.exists(traj_path):
        traj_path = os.mkdir(traj_path)

    mol_path = structure_dir + "/" + file_base + ".xyz"

    # check .xyz structure file
    #if not os.path.exists(mol_path):
    #    pdb2xyz(structure_dir, file_base)

    # Read molecule from xyz
    mol = read(mol_path)

    mol.set_calculator(ANIENS(aniensloader(model_ani2x_path, 0)))

    start_time = time.time()
    #dyn = LBFGS(mol, trajectory="%s/%s.traj" %(traj_path, file_base))
    dyn = LBFGS(mol)
    dyn.run(fmax=thr_fmax)
    run_time = time.time() - start_time
    print('[ANI Optimization - Total time:', run_time, 'seconds]')
    
    # Calculate energy
    ef = mol.get_potential_energy()
    print("Final Energy: ", ef)
    print("job is end for %s" %file_base)

    vib = Vibrations(mol)
    vib.run()
    return [ef, vib.get_energies()]

def termo_hormonic_limit(vib_energies, pot_energy):
    return HarmonicThermo(vib_analysis, pot_energy)

def runOptMultiMol(directory, base, seq_start, seq_end, thr_fmax):
    csv_file_name = "%s_opt_energies.csv"%directory
    list_results = []
    structure_dir = os.getcwd() + "/" + directory

    #file_names = [item for item in os.listdir(structure_dir) if ".pdb" in item]
    #file_names = [item for item in file_names if "A15B15" in item]
    #file_names = ["A14_B12.pdb"]
    calculated_files_path = structure_dir + "/" + csv_file_name
    if os.path.exists(calculated_files_path):
        df_calculated_files = pd.read_csv(calculated_files_path, index_col=0)
        calculated_files = df_calculated_files["Structure"].to_list()
    else:
        df_calculated_files = pd.DataFrame(list_results, columns=["Structure", "Optimzation_E_(eV)", "Run_time_(s)"])
        df_calculated_files.to_csv(calculated_files_path)
        calculated_files = df_calculated_files["Structure"].to_list()

    for i in range(seq_start, seq_end):
        #file_base = file_name.replace(".pdb", "")
        file_base = "{}{}".format(base, i)
        if file_base in calculated_files:
            #print(df_calculated_files.loc[df_calculated_files["Structure"] == file_base].iloc[0])
            list_results.append(df_calculated_files.loc[df_calculated_files["Structure"] == file_base].iloc[0].to_list())
            df_results = pd.DataFrame(list_results, columns=["Structure", "Optimzation_E_(eV)", "Run_time_(s)"])
            df_results.to_csv(calculated_files_path)
            print("Exist Energy for %s structure" %file_base)
            continue

        print("Calcultion starts for %s" %file_base)
        #convert from pdb to xyz
        #if not os.path.exists(mol_path):
        #    pdb2xyz(structure_dir, file_base)
        #try:
        results = runOptSingleMol(directory, file_base, thr_fmax)
        #except:
            #continue

        results.insert(0, file_base)

        list_results.append(results)
        df_results = pd.DataFrame(list_results, columns=["Structure", "Optimzation_E_(eV)", "Run_time_(s)"])
        df_results.to_csv("%s/%s" %(structure_dir, csv_file_name), mode="a", header=False)

#for i in range(5, 5):
#print("calculating ACS_M19_md_%s_solv"%i)


if args.calcMode == "sp_single_mol":
    #  file_base = "/M23_less_10001"
    structDIR = args.struct_dir
    file_base = args.fbase
    runSP(structDIR, file_base)
elif args.calcMode == "sp_multi_mol":
    structDIR = args.struct_dir
    namebase = args.namebase
    seq_start = args.seq_start
    seq_end = args.seq_end
    runSPmultiMol(structDIR, namebase, seq_start, seq_end)
elif args.calcMode == "sp_grouped_multi_mol":
    structDIR = args.struct_dir
    namebase = args.namebase
    seq_start = args.seq_start
    seq_end = args.seq_end
    group1= args.group1
    group2= args.group2
    runSPgroupedMultiMol(group1, group2, structDIR, namebase, seq_start, seq_end)
elif args.calcMode == "opt_single_mol":
    #  file_base = "/M23_less_10001"
    structDIR = args.struct_dir
    file_base = args.fbase
    thr_fmax= args.thr_fmax
    runOptSingleMol(structDIR, file_base, thr_fmax)

elif args.calcMode == "opt_multi_mol":
    #  file_base = "/M23_less_10001"
    structDIR = args.struct_dir
    namebase = args.namebase
    thr_fmax= args.thr_fmax
    seq_start = args.seq_start
    seq_end = args.seq_end
    runOptMultiMol(structDIR, namebase, seq_start, seq_end, thr_fmax)


elif args.calcMode == "opt_freq_single_mol":
    #  file_base = "/M23_less_10001"
    structDIR = args.struct_dir
    file_base = args.fbase
    thr_fmax= args.thr_fmax
    runOptFreqSingleMol(structDIR, file_base, thr_fmax)

