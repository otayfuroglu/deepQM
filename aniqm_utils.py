#
# Import ANI ensemble loader
import sys
#sys.path.append('/home/olexandr/notebooks/ASE_ANI/lib')
from ase_interface import ANIENS
from ase_interface import aniensloader

import  ase
import time
from ase import units
from ase.io import read, write, trajectory
from ase.optimize import BFGS, LBFGS

import pandas as pd
import os, re
from operator import itemgetter

import numpy  as np

def pdb2xyz(structure_path, file_base):

    fixed_pdb(structure_path, file_base)

    cmd = "obabel %s/%s.pdb -O %s/%s.xyz" %(structure_path, file_base, structure_path, file_base)
    os.system("%s" %cmd)
 
def fixed_atom_type(line):
    atom_type = "".join(re.findall("[a-zA-Z]+", line[11:17]))
    if atom_type.lower()[0:2] in ["cl", "br"]:
        line = line[:-6] + atom_type[:2] + " " + "\n"
    elif atom_type.lower()[0] in ["h", "c", "n", "o", "f", "i", "b","s"]:
        line = line[:-4] + atom_type[:1] + "  " + "\n"
    return line


def fixed_pdb(work_path, file_base, grp1, grp2):
    data_xyz = []
    lines = open("{}/{}.pdb".format(work_path, file_base)).readlines()
    fl_fixed_pdb = open("{}/{}.pdb".format(work_path, file_base), "w")
    for line in lines:
        if grp1 in line or grp2 in line:
            line=line.replace("1+","  ").replace("1-","  ")
            if line[-2:-1] == " " :
                line = fixed_atom_type(line)
        fl_fixed_pdb.write(line)
    fl_fixed_pdb.close()


def prepare_xyz_complex(work_path, file_base, grp1, grp2):

    fixed_pdb(work_path, file_base, grp1, grp2)

    fl_xyz = open("{}/{}.xyz".format(work_path, file_base), "w")
    #fl_xyz = open("{}.xyz".format(key_word), "w")
    data_xyz = []
    lines = open("{}/{}.pdb".format(work_path, file_base)).readlines()
    for line in lines:
        #print(line)
        if grp1 in line or grp2 in line:
            data_xyz.append(line)
    init_line = len(data_xyz)
    fl_xyz.write(str(init_line)+"\n{}/{}.xyz\n".format(work_path, file_base)) # file root to second row
    for line in data_xyz:
        split_line = line[30:].split()
        atom_sym = split_line[-1] # extract atom symbol 
        line = "\t".join(itemgetter( 0, 1, 2,)(split_line)) # select data for xyz format and convert from tuple to str
        line = atom_sym + "\t" + line


        fl_xyz.write(line)
        fl_xyz.write("\n")
    fl_xyz.close()

def prepare_xyz_grp(key_word, work_path, file_base):
    
    fl_xyz = open("{}/{}_{}.xyz".format(work_path, file_base, key_word), "w")
    #fl_xyz = open("{}.xyz".format(key_word), "w")
    data_xyz = []
    lines = open("{}/{}.pdb".format(work_path, file_base)).readlines()
    for line in lines:
        if key_word in line:
            data_xyz.append(line)
    init_line = len(data_xyz)
    fl_xyz.write(str(init_line)+"\n{}/{}_{}.xyz\n".format(work_path, file_base, key_word)) # file root to second row
    for line in data_xyz:
        split_line = line[30:].split()
        atom_sym = split_line[-1] # extract atom symbol 
        line = "\t".join(itemgetter( 0, 1, 2,)(split_line)) # select data for xyz format and convert from tuple to str
        line = atom_sym + "\t" + line

        fl_xyz.write(line)
        fl_xyz.write("\n")
    fl_xyz.close()
