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
parser.add_argument("-a", "--alpha", type=float, required=True,
                    help="Enter alpha constant")
parser.add_argument("-b", "--beta", type=float, required=True,
                    help="Enter beta constant")
args = parser.parse_args()


def avg(np_array):
    return np_array.mean()
    #  return df_mean


def std(np_array):
    return np_array.std()

if __name__ == "__main__":
    csv_path = args.csv_path
    alpha = args.alpha
    beta = args.beta
    column_names = ["diff_ani2x","diff_dftd3"]

    diff_ani2x = pd.read_csv(csv_path)["diff_ani2x"].to_numpy()
    diff_dftd3 = pd.read_csv(csv_path)["diff_dftd3"].to_numpy()

    diff_ani2x = diff_ani2x * EV2KT * KT2KCAL # kcal/mol
    diff_dftd3 = diff_dftd3 * EV2KT * KT2KCAL # kcal/mol

    diff_ani2xEn = avg(diff_ani2x)
    diff_dftd3En = avg(diff_dftd3)

    diff_ani2xStd = std(diff_ani2x)
    diff_dftd3Std = std(diff_dftd3)

    data_bindEn = diff_ani2x * alpha  + diff_dftd3 * beta
    bindEn = avg(data_bindEn)
    bindEnStd = std(data_bindEn)

    print(bindEn, bindEnStd)

    print("="*61)
    print("{:^61s}".format("SUMMARY"))
    print("="*61)
    print("DFTD3 (wb97x)            =   {:10.6f}    +/-     {:10.6f}".format(diff_ani2xEn, diff_ani2xStd))
    print("ANI-2x (wb97x/6-31G*)    =   {:10.6f}    +/-     {:10.6f}".format(diff_dftd3En, diff_dftd3Std))
    print("-"*61)
    print("Binding energy           =   {:10.6f}    +/-     {:10.6f}".format(bindEn, bindEnStd))
    print("-"*61)
    print("="*61)
