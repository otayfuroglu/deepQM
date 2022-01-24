#
import  ase
import time
from ase import units
from ase.io import read, write, trajectory

import pandas as pd
import os, re

import numpy  as np

# range of constants of index in pdb file
RECORDNAME = (0, 6)
ATOMTYPE = (12, 16)
ATOMSYM = (76, 78)
STARTCOORD = 30
LENCOORD = 8

def fixed_atom_type(line):
    atom_type = "".join(re.findall("[a-zA-Z]+", line[ATOMTYPE[0]:ATOMTYPE[1]]))
    if atom_type.lower()[0:2] in ["cl", "br", "na"]:
        line = line[:ATOMSYM[0]] + atom_type[:2] + " " + "\n"
    elif atom_type.lower()[0] in ["h", "c", "n", "o", "f", "i", "b","s"]:
        line = line[:ATOMSYM[0]] + atom_type[:1] + "  " + "\n"
    return line


def get_all_coords(structure_dir, file_base):
    coords = []
    lines = open("{}/{}.pdb".format(structure_dir, file_base)).readlines()
    for line in lines:
        if "ATOM" in line[RECORDNAME[0]:RECORDNAME[1]] or "HETATM" in line[RECORDNAME[0]:RECORDNAME[1]]:
            if line[ATOMSYM[0]:ATOMSYM[1]].strip() == "":
                # fiex missing atoms type end of column
                line = fixed_atom_type(line)
            coords.append(line)
    return coords

def coords2file(fl_xyz, coords):
    # coords start at 30th index in pdb file
    for line in coords:
        # atoms sym in 77-79 index in pdb file
        atom_sym = line[ATOMSYM[0]:ATOMSYM[1]].strip()
        coord = [line[STARTCOORD+pointer:STARTCOORD+pointer+LENCOORD].strip()
                 for pointer in [0, LENCOORD, 2*LENCOORD]] # for x, y, x coordinates
        line = "\t".join(coord) # select data for xyz format and convert from tuple to str
        line = atom_sym + "\t" + line

        fl_xyz.write(line)
        fl_xyz.write("\n")
    fl_xyz.close()

def prepare_xyz_files_grouped(structure_dir, file_base, index_file_path, grp1, grp2):

    sys_coords = []
    for grps in [grp1, grp2]:
        grp_coords = []
        for grp in set(grps):
            grp_coords += get_grp_coords(grp, structure_dir, index_file_path, file_base).tolist()
        sys_coords += grp_coords

        grp_xyz_path = "{}/{}_grp_{}.xyz".format(structure_dir, file_base, "".join(map(str, grps)))
        fl_xyz_grp = open(grp_xyz_path, "w")

        init_line = len(grp_coords)

        fl_xyz_grp.write(str(init_line)+"\n\n")
        coords2file(fl_xyz_grp, grp_coords)

    fl_xyz_sys = open("{}/{}.xyz".format(structure_dir, file_base), "w")
    init_line = len(sys_coords)
    fl_xyz_sys.write(str(init_line)+"\n{}/{}.xyz\n".format(structure_dir, file_base)) # file root to second row
    coords2file(fl_xyz_sys, sys_coords)


def prepare_xyz_files(structure_dir, file_base):

    sys_coords = get_all_coords(structure_dir, file_base)
    fl_xyz_sys = open("{}/{}.xyz".format(structure_dir, file_base), "w")
    init_line = len(sys_coords)
    fl_xyz_sys.write(str(init_line)+"\n{}/{}.xyz\n".format(structure_dir, file_base)) # file root to second row
    coords2file(fl_xyz_sys, sys_coords)


def get_grp_coords(grp, structure_dir, index_file_path, file_base):
    lines = open(index_file_path,"r").readlines()

    index_groups = []
    for i, line in enumerate(lines):
        if "[ " in line:
            index_groups += [i]

    assert grp + 1 <= len(index_groups), "Out of system index size"
    index_group = index_groups[grp:grp+2]

    a = index_group[0]+1
    try:
        grp_lines = lines[a:index_group[1]]
    except:
        grp_lines = lines[a:]

    grp_index = []
    for items in grp_lines:
        grp_index += [int(item) for item in items.split()]

    grp_index = np.array(grp_index) - 1
    coords = np.array(get_all_coords(structure_dir, file_base))
    return coords[grp_index]

