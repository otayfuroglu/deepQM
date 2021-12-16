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
from aniqm_utils import prepare_xyz_grp, prepare_xyz_complex, pdb2xyz
import torchani

import torch
print("Nuber of CUDA devices: ", torch.cuda.device_count())
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-aniDIR", "--aniDIR",
                    type=str, required=True,
                    help="")
parser.add_argument("-calcMode", "--calcMode",
                    type=str, required=True,
                    help="")
parser.add_argument("-model_list", "--model_list", nargs='+', default=[],
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
        labels +=  ["SPE_%s_%s_%s"%(model_name, grp1, grp2),
                   "SPE_%s_%s"%(model_name, grp1),
                   "SPE_%s_%s"%(model_name, grp2), "diff_%s" %model_name]

    df_results = pd.DataFrame(columns=labels)
    df_results.to_csv("%s/%s_SP_energies.csv" %(structure_dir, base))

    for i in range(seq_start, seq_end):
        file_base = "{}{}".format(base, i)

        prepare_xyz_complex(structure_dir, file_base, grp1, grp2)
        mol_path = structure_dir + "/" + file_base + ".xyz"
        mol = read(mol_path)

        data = {model_name: [] for model_name in model_names}
        print("%d. pdb file is processing ..." %i)
        #file_base = file_name.replace(".pdb", "").replace(".xyz", "")
        print("Calcultion starts for %s" %file_base)

        j = 0
        for model_name, model in load_models(model_names).items():
            result = calcSPWithModel(model, mol)
            data[f"{model_name}"].append(calcSPWithModel(model, mol))
            if result == 0.0:
                for grp in [grp1, grp2]:
                    data[f"{model_name}"].append(result)
                continue
            else:
                for grp in [grp1, grp2]:
                    if j == 0:
                        prepare_xyz_grp(grp, structure_dir, file_base)
                        file_base_new = "{}_{}".format(file_base, grp)
                        mol_path = structure_dir + "/" + file_base_new + ".xyz"
                        mol = read(mol_path)
                    data[f"{model_name}"].append(calcSPWithModel(model, mol))
            j += 1

        list_results = [file_base]
        for model_name in model_names:
            model_data = data[model_name]
            list_results += model_data
            list_results += get_diff(model_data)
        df_results = pd.DataFrame([list_results], columns=labels)
        df_results.to_csv("%s/%s_SP_energies.csv" %(structure_dir, base), mode="a", header=False)


args = parser.parse_args()
model_names = args.model_list

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

