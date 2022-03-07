#
import  ase
import time
from ase import units
from ase.io import read, write, trajectory
from ase.optimize import BFGS, LBFGS
from ase.vibrations import Vibrations

import pandas as pd
import os, re, sys
from operator import itemgetter
import tqdm

import shutil
import numpy  as np
import argparse
#  from multiprocessing import Pool
from torch.multiprocessing import Pool, Process, set_start_method
from multiprocessing import current_process

#  from aimnet import load_AIMNetMT_ens, load_AIMNetSMD_ens, AIMNetCalculator
from deepqm_utils import prepare_xyz_files_grouped, prepare_xyz_files
import torchani

import torch
ngpu = torch.cuda.device_count()
print("Nuber of CUDA devices: ", ngpu)


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("calcMode", type=str)
parser.add_argument("n_procs", type=int)
parser.add_argument("model_list", nargs='+', default=[], type=str)
parser.add_argument("struct_dir", type=str)
parser.add_argument("namebase", type=str, )
parser.add_argument("seq_start", type=int, default=1)
parser.add_argument("seq_end", type=int, default=1E6)
parser.add_argument("index_file_path", type=str)
parser.add_argument("group1", type=str)
parser.add_argument("group2", type=str)
parser.add_argument("thr_fmax", type=float, default=0.01)
parser.add_argument("maxiter", type=float, default=500)



warning = """Warning: an error was encountered.
CUDA out of memory in case of aimnet model or your system may contain elements such as \"F, Cl and S\".
Returned 0.0 for energy"""

def calcSPWithModel(calculator, mol):

    mol.set_calculator(calculator)

    # check atom type for ani2x
    #  atom_types = mol.get_chemical_symbols()
    #  ani2_atom_type = [sym for sym in atom_types if sym in ["Cl", "F", "S"]]
    #  if len(ani2_atom_type) >= 1:
    #      return 0.0

    # try:
    return np.round(mol.get_potential_energy(), 6)
    # except:
    #     print(warning)
    #     return 0.0


def getD3calc(xc="pbe"):
    from ase.calculators.dftd3 import DFTD3
    procid = current_process().pid
    workdir = f"ase_dft{procid}"

    # due to ase bug, chdir working dir and comment out label and directory
    pwd = os.getcwd()
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.mkdir(workdir)
    os.chdir(workdir)

    params = {'s6': 1.0, 's9': 1.0, 'alp': 16.0, 's8': 0.6633, 'a1': 0.4288, 'a2': 3.9935}
    return DFTD3(
        #  label="tmp_d3",
        #  directory=procid,
        xc=xc
        #  damping="bj",
        #  a1=0.4288,
        #  a2=3.9935,
        #  s6=1.0,
        #  #  s9=1.0,
        #  alpha6=16.0,
        #  s8=0.6633,
        
    )
    os.chdir(pwd)
    shutil.rmtree("ase_dft*") #intesting run after retunr statement


def getD4calc(xc="pbe"):
    from dftd4.ase import DFTD4
    procid = current_process()
    return DFTD4(
        label="tmp_d3_%s" % procid.pid,
        method=xc
    )


def setG16Calculator():
    from ase.calculators.gaussian import Gaussian
    procid = current_process()

    calculator = Gaussian(
        label="tmp_g16/%s" % procid.pid,
        #  chk="tmp.chk",
        xc="wb97x",
        basis="6-31g*",
        scf="maxcycle=100",
    )
    return calculator


def load_calculators(model_names, device):
    # model_list ani
    models = {}
    model_names = [model_name.lower() for model_name in model_names]
    if "ani1x" in model_names:
        models["ani1x"] = torchani.models.ANI1x().to(device).ase()
    if "ani1ccx" in model_names:
        models["ani1ccx"] = torchani.models.ANI1ccx().to(device).ase()
    if "ani2x" in model_names:
        models["ani2x"] = torchani.models.ANI2x().to(device).ase()
    if "dftd3" in model_names:
        if shutil.which("dftd3") is None:
            print("dftd3 program CAN NOT FOUND !!!")
            print("You can install dftd3 by running the following command in command line")
            print("conda install -c psi4 dftd3")
            sys.exit(1)
        else:
            models["dftd3"] = getD3calc()
    if "dftd4" in model_names:
        if shutil.which("dftd4") is None:
            print("dftd4 program CAN NOT FOUND !!!")
            print("You can install dftd4 by running the following command in command line")
            print("conda install -c conda-forge dftd4")
            sys.exit(1)
        else:
            models["dftd4"] = getD4calc()
    if "g16" in model_names:
        if shutil.which("g16") is None:
            print("g16 program CAN NOT FOUND !!!")
            sys.exit(1)
        else:
            models["g16"] = setG16Calculator()
    #  if "dftd4" in model_names:

    #      models["dftd4"] = getD4calc()
    # else:
    #     print("All calcultors or one CAN NOT FOUND !!!")
    #     sys.exit(1)

    # model list aimnet
    #  model_gas = load_AIMNetMT_ens().to(device)
    #  model_smd = load_AIMNetSMD_ens().to(device)
    #  aimnet_gas = AIMNetCalculator(model_gas)
    #  aimnet_smd = AIMNetCalculator(model_smd)

    #  models = {"ani1x": ani1x, "ani1ccx": ani1ccx, "ani2x": ani2x,}
              #  "aimnetgas": aimnet_gas, "aimnetsmd": aimnet_smd}

    #  return {key: models[key] for key in model_names}
    return models


def _getLabels():

    labels = ["Structure"]
    for model_name in model_names:
        labels +=  ["%sE_%s"%(calc_type, model_name)]
    return labels


def _getLabelsGrouped():

    labels = ["Structure"]
    for model_name in model_names:
        labels +=  ["%sE_%s_%s_%s"%(calc_type, model_name, "_".join(map(str, grp1)), "_".join(map(str, grp2))),
                   "%sE_%s_%s"%(calc_type, model_name, "_".join(map(str, grp1))),
                   "%sE_%s_%s"%(calc_type, model_name, "_".join(map(str, grp2))), "diff_%s" %model_name]
    return labels


def _getcsvPathGrouped() :
    return "%s/%s_%s_energies_%s_%s.csv" %(structure_dir, namebase, calc_type,
                                "_".join(map(str, grp1)),
                                "_".join(map(str, grp2)))


def _getcsvPath() :
    return "%s/%s_%s_energies.csv" %(structure_dir, namebase, calc_type)


def _getFilenames():
    if seq_end == -1:
        return [item for item in os.listdir(structure_dir) if ".pdb" in item]
    else:
        return ["{}{}.pdb".format(namebase, i) for i in range(seq_start, seq_end)]


def _getDevices(idx):
    if ngpu != 0:
        #  os.environ["CUDA_VISIBLE_DEVICES"] = str(idx % ngpu)
        return torch.device('cuda:%s' %str(idx % ngpu))
    else:
        return torch.device("cpu")


def get_diff(model_data):
    return [model_data[0] - model_data[1] - model_data[2]]


def checkZeroError(data):
    if sum([sum(list_val) for list_val in data.values()]) == 0.0:
        print("\n!!! All calculation result returned 0.0 !!!")
        print("May be CUDA out of memory error.. Consider number of atoms in system")
        return True
    else:
        return False


def _SPgroupedMultiMol(idx):
    device = _getDevices(idx)

    file_names = _getFilenames()
    file_name = file_names[idx]
    file_base = file_name.replace(".pdb", "")

    prepare_xyz_files_grouped(structure_dir, file_base, index_file_path, grp1, grp2)
    mol_path = structure_dir + "/" + file_base + ".xyz"
    mol = read(mol_path)

    data = {model_name: [] for model_name in model_names}
    print("%s. pdb file is processing ..." %file_name)

    for model_name, model in load_calculators(model_names, device).items():
        result = calcSPWithModel(model, mol)
        data[f"{model_name}"].append(result)
        if result == 0.0:
            for grp in [grp1, grp2]:
                data[f"{model_name}"].append(result)
            continue
        else:
            for grp in [grp1, grp2]:
                #  for grp in set(grps):
                file_base_new = "{}_grp_{}".format(file_base, "".join(map(str, grp)))
                mol_path = structure_dir + "/" + file_base_new + ".xyz"
                mol_grp = read(mol_path)
                result = calcSPWithModel(model, mol_grp)
                data[f"{model_name}"].append(result)

    list_results = [file_base]
    for model_name in model_names:
        model_data = data[model_name]
        list_results += model_data
        list_results += get_diff(model_data)
    #  print(set(grp1))
    labels = _getLabelsGrouped()
    csv_path = _getcsvPathGrouped()
    df_results = pd.DataFrame([list_results], columns=labels)
    df_results.to_csv(csv_path, mode="a", header=False)

    # for check sum of Zero error
    return data


def runSPgroupedMultiMol(n_procs):
    labels = ["Structure"]

    labels = _getLabelsGrouped()
    df_results = pd.DataFrame(columns=labels)

    csv_path = _getcsvPathGrouped()
    df_results.to_csv(csv_path)

    file_names = _getFilenames()
    idxs = range(len(file_names))

    # run on multiprocessors
    pool = Pool(processes=n_procs)
    result_list_tqdm = []

    # implementation of  multiprocessor in tqdm. Ref.https://leimao.github.io/blog/Python-tqdm-Multiprocessing/
    for result in tqdm.tqdm(pool.imap_unordered(func=_SPgroupedMultiMol, iterable=idxs), total=len(idxs)):
        #  check sum of Zero error
        if checkZeroError(result):
            sys.exit(1)
        result_list_tqdm.append(result)


def _SPMultiMol(idx):
    device = _getDevices(idx)

    file_names = _getFilenames()
    file_name = file_names[idx]
    file_base = file_name.replace(".pdb", "")

    prepare_xyz_files(structure_dir, file_base)
    mol_path = structure_dir + "/" + file_base + ".xyz"
    mol = read(mol_path)

    data = {model_name: [] for model_name in model_names}
    print("%s. pdb file is processing ..." %file_name)

    for model_name, model in load_calculators(model_names, device).items():
        result = calcSPWithModel(model, mol)
        data[f"{model_name}"].append(result)

    list_results = [file_base]
    for model_name in model_names:
        model_data = data[model_name]
        list_results += model_data
    #  print(set(grp1))
    labels = _getLabels()
    csv_path = _getcsvPath()
    df_results = pd.DataFrame([list_results], columns=labels)
    df_results.to_csv(csv_path, mode="a", header=False)

    # for check sum of Zero error
    return data


def runSPMultiMol(n_procs):
    labels = ["Structure"]

    labels = _getLabels()
    df_results = pd.DataFrame(columns=labels)

    csv_path = _getcsvPath()
    df_results.to_csv(csv_path)

    file_names = _getFilenames()
    idxs = range(len(file_names))

    # run on multiprocessors
    pool = Pool(processes=n_procs)
    result_list_tqdm = []

    # implementation of  multiprocessor in tqdm. Ref.https://leimao.github.io/blog/Python-tqdm-Multiprocessing/
    for result in tqdm.tqdm(pool.imap_unordered(func=_SPMultiMol, iterable=idxs), total=len(idxs)):
        #  check sum of Zero error
        if checkZeroError(result):
            sys.exit(1)
        result_list_tqdm.append(result)


def runOptSingleMol(file_base, device):

    mol_path = structure_dir + "/" + file_base + ".xyz"
    mol = read(mol_path)

    ani2x = torchani.models.ANI2x().to(device).ase()
    mol.set_calculator(ani2x)

    #  try:
    dyn = LBFGS(mol, logfile="opt.log")
    dyn.run(fmax=thr_fmax, steps=args.maxiter)

    # Calculate energy
    ef = mol.get_potential_energy()
    write("%s/optimized_%s.xyz" %(structure_dir, file_base), mol)
    return ef
    #  except:
    #      print(warning)
    #      return 0.0


def _optGroupedMultiMol(idx):
    device = _getDevices(idx)

    file_names = _getFilenames()
    file_name = file_names[idx]
    file_base = file_name.replace(".pdb", "")

    prepare_xyz_files_grouped(structure_dir, file_base, index_file_path, grp1, grp2)

    runOptSingleMol(file_base, device)
    opt_xyz_path = "%s/optimized_%s.xyz" %(structure_dir, file_base)
    opt_pdb_path = "%s/%s.pdb" %(structure_dir, file_base)
    #  os.remove(opt_pdb_path)
    os.system("obabel %s -O %s" %(opt_xyz_path, opt_pdb_path))


def runOptGroupedMultiMol(n_procs, thr_fmax):

    file_names = _getFilenames()
    idxs = range(len(file_names))

    # to optimizing strucututes
    # for waiting all multi processes ends
    results= []
    with Pool(n_procs) as pool:
        result = pool.map_async(_optGroupedMultiMol, idxs)
        pool.close()
        pool.join()
        results.append(result)
    # wait for done all processes
    [result.wait() for result in results]

    # to run sp grouped
    runSPgroupedMultiMol(n_procs)


def _OptMultiMol(idx):

    device = _getDevices(idx)

    file_names = _getFilenames()
    file_name = file_names[idx]
    file_base = file_name.replace(".pdb", "")

    data = {model_name: [] for model_name in model_names}
    print("%s. pdb file is processing ..." %file_name)

    for model_name, model in load_calculators(model_names, device).items():
        result = runOptSingleMol(file_base, device)
        data[f"{model_name}"].append(result)

    list_results = [file_base]
    for model_name in model_names:
        model_data = data[model_name]
        list_results += model_data
    #  print(set(grp1))
    labels = _getLabels()
    csv_path = _getcsvPath()
    df_results = pd.DataFrame([list_results], columns=labels)
    df_results.to_csv(csv_path, mode="a", header=False)

    # for check sum of Zero error
    return data


def runOptMultiMol(n_procs, thr_fmax):

    file_names = _getFilenames()
    idxs = range(len(file_names))

    # initialize csv file
    labels = _getLabels()
    df_results = pd.DataFrame(columns=labels)

    csv_path = _getcsvPath()
    df_results.to_csv(csv_path)

    #  # to optimizing strucututes
    #  with Pool(n_procs) as pool:
    #      pool.map(_OptMultiMol, idxs)
    # implementation of  multiprocessor in tqdm. Ref.https://leimao.github.io/blog/Python-tqdm-Multiprocessing/
    for result in tqdm.tqdm(pool.imap_unordered(func=_OptMultiMol, iterable=idxs), total=len(idxs)):
        #  check sum of Zero error
        if checkZeroError(result):
            sys.exit(1)
        result_list_tqdm.append(result)


args = parser.parse_args()
n_procs = args.n_procs
model_names = args.model_list
structure_dir = args.struct_dir
namebase = args.namebase
seq_start = args.seq_start
seq_end = args.seq_end
grp1 =  [int(i) for i in args.group1.split("_")]
grp2 =  [int(i) for i in args.group2.split("_")]

#rm .xyz file in structure_dir

if "grouped" in args.calcMode:
    index_file_path = args.index_file_path
    if not os.path.isfile(index_file_path):
        print("Error: Not found index file")
        sys.exit(1)

if "sp" in args.calcMode.lower():
    calc_type = "SP"
else:
    calc_type = "Opt"
    thr_fmax= args.thr_fmax

if __name__ == "__main__":
    try:
         set_start_method('spawn')
    except RuntimeError:
        pass

    if args.calcMode == "sp_single_mol":
        n_procs = 1
        runSPMultiMol(n_procs)
    elif args.calcMode == "sp_multi_mol":
        runSPMultiMol(n_procs)
    elif args.calcMode == "sp_grouped_multi_mol":
        runSPgroupedMultiMol(n_procs)
    elif args.calcMode == "opt_grouped_multi_mol":
        runOptGroupedMultiMol(n_procs, thr_fmax)
    elif args.calcMode == "opt_single_mol":
        n_procs = 1
        runOptMultiMol(n_procs, thr_fmax)
    elif args.calcMode == "opt_multi_mol":
        runOptMultiMol(n_procs, thr_fmax)
