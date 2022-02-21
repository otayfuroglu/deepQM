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
                    help="Enter alpha constant")
parser.add_argument("-b", "--beta", type=float, default=0.0, required=False,
                    help="Enter beta constant")
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
    if beta != 0.0:
        if "diff_dftd3" not in column_names:
            print("Error: Not found dftd3 in methods")
            sys.exit()

        diff_dftd3 = df["diff_dftd3"].to_numpy()
        diff_dftd3 = diff_dftd3 * EV2KT * KT2KCAL # kcal/mol
        diff_dftd3En = avg(diff_dftd3)
        diff_dftd3Std = std(diff_dftd3)

    data_bindEn = diff_ani * alpha  + diff_dftd3 * beta
    bindEn = avg(data_bindEn)
    bindEnStd = std(data_bindEn)

    print("="*61)
    print("{:^61s}".format("SUMMARY"))
    print("="*61)
    if beta != 0.0:
        print("DFTD3 (wb97x)            =   {:10.6f}    +/-     {:10.6f}".format(diff_dftd3En, diff_dftd3Std))
    if alpha != 1.0:
        print("{} (wb97x/6-31G*)    =   {:10.6f}    +/-     {:10.6f}".format(ani_method, diff_aniEn, diff_aniStd))
    print("-"*61)
    print("Binding energy           =   {:10.6f}    +/-     {:10.6f}".format(bindEn, bindEnStd))
    print("-"*61)
    print("="*61)
