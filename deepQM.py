#
import sys

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

#  from aimnet import load_AIMNetMT_ens, load_AIMNetSMD_ens, AIMNetCalculator
from deepqm_utils import prepare_xyz_files
import torchani

import torch
print("Nuber of CUDA devices: ", torch.cuda.device_count())
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("calcMode", type=str)
parser.add_argument("model_list", nargs='+', default=[], type=str)
parser.add_argument("struct_dir", type=str)
parser.add_argument("fbase", type=str)
parser.add_argument("namebase", type=str, )
parser.add_argument("seq_start", type=int, default=1)
parser.add_argument("seq_end", type=int, default=1E6)
parser.add_argument("index_file_path", type=str)
parser.add_argument("group1", type=str)
parser.add_argument("group2", type=str)
parser.add_argument("thr_fmax", type=float, default=0.01)


def calcSPWithModel(calculator, mol):

    mol.set_calculator(calculator)

    # check atom type for ani2x
    #  atom_types = mol.get_chemical_symbols()
    #  ani2_atom_type = [sym for sym in atom_types if sym in ["Cl", "F", "S"]]
    #  if len(ani2_atom_type) >= 1:
    #      return 0.0

    worning = """Warning: an error was encountered.
    CUDA out of memory in case of aimnet model or your system may contain elements such as \"F, Cl and S\".
    Returned 0.0 for energy"""
    try:
        return np.round(mol.get_potential_energy(), 8)
    except:
        print(worning)
        return 0.0


def load_models(model_names):
    # model_list ani
    ani1x = torchani.models.ANI1x().to(device).ase()
    ani1ccx = torchani.models.ANI1ccx().to(device).ase()
    ani2x = torchani.models.ANI2x().to(device).ase()

    # model list aimnet
    #  model_gas = load_AIMNetMT_ens().to(device)
    #  model_smd = load_AIMNetSMD_ens().to(device)
    #  aimnet_gas = AIMNetCalculator(model_gas)
    #  aimnet_smd = AIMNetCalculator(model_smd)

    models = {"ani1x": ani1x, "ani1ccx": ani1ccx, "ani2x": ani2x,}
              #  "aimnetgas": aimnet_gas, "aimnetsmd": aimnet_smd}

    return {key: models[key] for key in model_names}


def get_diff(model_data):
    return [model_data[0] - model_data[1] - model_data[2]]


def runSPgroupedMultiMol(grp1, grp2, directory, base, seq_start, seq_end):
    structure_dir = directory
    labels = ["Structure"]

    for model_name in model_names:
        labels +=  ["%sE_%s_%s_%s"%(calc_type, model_name, grp1, grp2),
                   "%sE_%s_%s"%(calc_type, model_name, grp1),
                   "%sE_%s_%s"%(calc_type, model_name, grp2), "diff_%s" %model_name]

    df_results = pd.DataFrame(columns=labels)

    csv_path = "%s/%s_%s_energies_%s_%s.csv" %(structure_dir, base, calc_type,
                                "_".join(map(str, grp1)),
                                "_".join(map(str, grp2)))
    df_results.to_csv(csv_path)

    if seq_end == 0:
        file_names = [item for item in os.listdir(structure_dir) if ".pdb" in item]
    else:
        file_names = ["{}{}.pdb".format(base, i) for i in range(seq_start, seq_end)]

    #  for i in range(seq_start, seq_end):
    for file_name in file_names:
        file_base = file_name.replace(".pdb", "")

        prepare_xyz_files(structure_dir, file_base, args.index_file_path, grp1, grp2)
        mol_path = structure_dir + "/" + file_base + ".xyz"
        mol = read(mol_path)

        data = {model_name: [] for model_name in model_names}
        print("%s. pdb file is processing ..." %file_name)

        j = 0
        for model_name, model in load_models(model_names).items():
            result = calcSPWithModel(model, mol)
            data[f"{model_name}"].append(calcSPWithModel(model, mol))
            if result == 0.0:
                for grp in [grp2, grp2]:
                    data[f"{model_name}"].append(result)
                continue
            else:
                for grp in [grp1, grp2]:
                    #  for grp in set(grps):
                    if j == 0:
                        #  prepare_xyz_grp(grp, structure_dir, args.index_file_path, file_base)
                        file_base_new = "{}_grp_{}".format(file_base, "".join(map(str, grp)))
                        mol_path = structure_dir + "/" + file_base_new + ".xyz"
                        mol = read(mol_path)
                    data[f"{model_name}"].append(calcSPWithModel(model, mol))
            j += 1

        list_results = [file_base]
        for model_name in model_names:
            model_data = data[model_name]
            list_results += model_data
            list_results += get_diff(model_data)
        print(set(grp1))
        df_results = pd.DataFrame([list_results], columns=labels)
        df_results.to_csv(csv_path, mode="a", header=False)


def runOptSingleMol(structure_dir, file_base, thr_fmax):

    mol_path = structure_dir + "/" + file_base + ".xyz"
    mol = read(mol_path)

    ani2x = torchani.models.ANI2x().to(device).ase()
    mol.set_calculator(ani2x)

    start_time = time.time()
    #dyn = LBFGS(mol, trajectory="%s/%s.traj" %(traj_path, file_base))
    dyn = LBFGS(mol)
    dyn.run(fmax=thr_fmax)
    run_time = time.time() - start_time
    # print('[ANI Optimization - Total time:', run_time, 'seconds]')

    # Calculate energy
    ef = mol.get_potential_energy()
    write("%s/optimized_%s.xyz" %(structure_dir, file_base), mol)
    #  print("Final Energy: ", ef)
    #  print("job is end for %s" %file_base)
    #  return [ef, run_time]
    #  return ef, mol


def runOptGroupedMultiMol(grp1, grp2, directory, base, seq_start, seq_end, thr_fmax):
    structure_dir = directory

    if seq_end == 0:
        file_names = [item for item in os.listdir(structure_dir) if ".pdb" in item]
    else:
        file_names = ["{}{}.pdb".format(base, i) for i in range(seq_start, seq_end)]

    #  for i in range(seq_start, seq_end):
    for file_name in file_names:
        file_base = file_name.replace(".pdb", "")
        prepare_xyz_files(structure_dir, file_base, args.index_file_path, grp1, grp2)

        runOptSingleMol(structure_dir, file_base, thr_fmax)
        opt_xyz_path = "%s/optimized_%s.xyz" %(structure_dir, file_base)
        opt_pdb_path = "%s/%s.pdb" %(structure_dir, file_base)
        #  os.remove(opt_pdb_path)
        os.system("obabel %s -O %s" %(opt_xyz_path, opt_pdb_path))

    runSPgroupedMultiMol(grp1, grp2, directory, base, seq_start, seq_end)


args = parser.parse_args()
model_names = args.model_list
print(model_names)

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
    calc_type = "SP"
    structDIR = args.struct_dir
    namebase = args.namebase
    seq_start = args.seq_start
    seq_end = args.seq_end
    group1 = [int(i) for i in args.group1.split("_")]
    group2 = [int(i) for i in args.group2.split("_")]
    runSPgroupedMultiMol(group1, group2, structDIR, namebase, seq_start, seq_end)
elif args.calcMode == "opt_grouped_multi_mol":
    #  file_base = "/M23_less_10001"
    calc_type = "Opt"
    directory = args.struct_dir
    base = args.namebase
    thr_fmax= args.thr_fmax
    seq_start = args.seq_start
    seq_end = args.seq_end
    grp1 = [int(i) for i in args.group1.split("_")]
    grp2 = [int(i) for i in args.group2.split("_")]
    runOptGroupedMultiMol(grp1, grp2, directory, base, seq_start, seq_end, thr_fmax)


elif args.calcMode == "opt_single_mol":
    #  file_base = "/M23_less_10001"
    structDIR = args.struct_dir
    file_base = args.fbase
    thr_fmax= args.thr_fmax
    runOptSingleMol(structDIR, file_base, thr_fmax)

elif args.calcMode == "opt_multi_mol":
    #  file_base = "/M23_less_10001"
    calc_type = "Opt"
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

