#
import os, sys
import pandas as pd
import numpy as np
import argparse


OUT_DIR = "./"
EV2KT = 38.94
KT2KCAL = 0.593

parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-in", "--csv_path", type=str, required=True,
                    help="Enter generated .csv file by deepQM")
parser.add_argument("-ani", "--ani_method", type=str, required=True,
                    help="Enter one of the ani versions (ani1x, ani1ccs or ani2x)")
parser.add_argument("-a", "--alpha", type=float, default=1.0, required=False,
                    help="Enter alpha coefficient for ani")
parser.add_argument("-b", "--beta", type=float, default=0.127, required=False,
                    help="Enter beta coefficient for d3")
parser.add_argument("-g", "--gamma", type=float, default=-5.111, required=False,
                    help="Enter gamma coefficient for intercept")
parser.add_argument("-ls_in_ls", "--ls_in_ls", type=float, default=0, required=False,
                    help="Enter gamma coefficient for intercept")
parser.add_argument("-ls_in_pls", "--ls_in_pls", type=float, default=0, required=False,
                    help="Enter gamma coefficient for intercept")
args = parser.parse_args()


def avg(np_array):
    return np_array.mean()
    #  return df_mean


def std(np_array):
    return np_array.std()

if __name__ == "__main__":
    csv_path = args.csv_path
    ani_method = args.ani_method
    alpha = args.alpha
    beta = args.beta
    gamma=args.gamma
    #sim_type = args.sim_type
    ls_in_ls = args.ls_in_ls
    ls_in_pls = args.ls_in_pls

    df = pd.read_csv(csv_path)
    column_names = df.columns

    if "diff_%s"%ani_method not in column_names:
        print("Error: Not found %s in ani methods" %ani_method)
        sys.exit()

    diff_ani = df["diff_%s" %ani_method].to_numpy()
    diff_ani = diff_ani * EV2KT * KT2KCAL # kcal/mol
    diff_aniEn = avg(diff_ani)
    diff_aniStd = std(diff_ani)


    diff_dftd3 = 0.0
    if alpha != 0.0:
        if "diff_dftd3" not in column_names:
            print("Error: Not found dftd3 in methods")
            sys.exit()

        diff_dftd3 = df["diff_dftd3"].to_numpy()
        diff_dftd3 = diff_dftd3 * EV2KT * KT2KCAL # kcal/mol
        diff_dftd3En = avg(diff_dftd3)
        diff_dftd3Std = std(diff_dftd3)


    if ls_in_ls == 0 or ls_in_pls == 0:
        data_bindEn = diff_dftd3 * alpha  + diff_ani * beta + gamma
    else:
        data_bindEn = (diff_dftd3 + ls_in_pls - ls_in_ls) * alpha  + (diff_ani + ls_in_pls - ls_in_ls) * beta + gamma

    bindEn = avg(data_bindEn)
    bindEnStd = std(data_bindEn)
    with open('Summary.dat','w') as f:

        print("="*61,file=f)
        print("{:^61s}".format("SUMMARY-kcal/mol"),file=f)
        print("="*61,file=f)
        if alpha != 0.0:
            print("dftd3 (wb97x)            =   {:10.6f}    +/-     {:10.6f}".format(diff_dftd3En, diff_dftd3Std),file=f)
        if beta != 0.127:
            print("{} (wb97x/6-31G*)    =   {:10.6f}    +/-     {:10.6f}".format(ani_method, diff_aniEn, diff_aniStd),file=f)
        if ls_in_ls != 0 and ls_in_pls != 0:
            print("ligand solvation term  from LS sim  =   {:10.6f}".format(ls_in_ls),file=f)
            print("ligand solvation term from PLS sim   =   {:10.6f}".format(ls_in_pls),file=f)
        print("-"*61,file=f)
        print("Binding energy           =   {:10.6f}    +/-     {:10.6f}".format(bindEn, bindEnStd),file=f)
        print("-"*61,file=f)
        print("="*61,file=f)
