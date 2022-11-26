#
import os
import tqdm
import argparse

arg_handler = argparse.ArgumentParser(description="")
arg_handler.add_argument("-cif_dir",  type=str, required=True)


def fixPymolcif(fl_path, fl_fixed_path):
    fl_fixed = open(fl_fixed_path, "w")
    with open(fl_path) as lines:
        for line in lines:
            if line.startswith("_"):
                sp_line = line.split(" ")
                fixed_part = sp_line[0].replace(".", "_")
                sp_line[0] = fixed_part
                line = " ".join(sp_line)
            fl_fixed.write(line)


if __name__ == "__main__":
    args = arg_handler.parse_args()
    cif_dir = args.cif_dir
    fl_fixed_dir = cif_dir.replace("/", "") + "_fixed"
    if not os.path.exists(fl_fixed_dir):
        os.mkdir(fl_fixed_dir)
    for fl_name in tqdm.tqdm([name for name in os.listdir(cif_dir)
                              if name.endswith(".cif")]):
        fl_path = os.path.join(cif_dir, fl_name)
        fl_fixed_path = os.path.join(fl_fixed_dir, fl_name)
        fixPymolcif(fl_path, fl_fixed_path)
