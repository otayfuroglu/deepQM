#
import os, sys
import pandas as pd
import numpy as np
import argparse

OUTPUT_DIR = "./"
EV2KT = 38.94
KT2KCAL = 0.593
method_list = ["avg", "exp", "jar", "cumu"]


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-OUTPUT_DIR", "--OUTPUT_DIR", type=str, required=True,
                    help="Enter output directory which has generated .csv file by deepQM")
args = parser.parse_args()


def avg(np_array):
    return np_array.mean()
    #  return df_mean

def exp(np_array):
    values = np_array * -1
    values = np.exp(values)
    values_mean = values.mean()
    result = np.log(values_mean) * -1
    return result

def jar(np_array):
    values_mean = np_array.mean()
    deviation = (np_array - values_mean) * -1
    exp_dev = np.exp(deviation)
    mean = exp_dev.mean()
    result = values_mean - np.log(mean) * -1
    return result

def cumu(np_array):
    values_mean = np_array.mean()
    values_var = np.var(np_array)
    result = values_mean - 0.5 * values_var
    return result


if __name__ == "__main__":
    OUTPUT_DIR = args.OUTPUT_DIR
    file_bases =[item.replace(".csv", "") for item in os.listdir(OUTPUT_DIR) if ".csv" in item]

    column_names = ["diff_ani2x","diff_dftd3"]
    df_list = []
    for column_name in column_names:
        for method in method_list:
            results = []
            for file_base in file_bases:
                csv_path = "%s/%s.csv" % (OUTPUT_DIR, file_base)
                df = pd.read_csv(csv_path)
                values = df[column_name].to_list() # in eV
                values = np.array(values) * EV2KT # in kT
                #  mean = values.mean()

                #  print(method)
                if method == "avg":
                    result_column_name = method + "_" + column_name
                    result = avg(values)
                    result = result * KT2KCAL # kcal/mol

                elif method == "exp":
                    result_column_name = method + "_" + column_name
                    result = exp(values)
                    result = result * KT2KCAL # kcal/mol

                elif method == "jar":
                    result_column_name = method + "_" + column_name
                    result = jar(values)
                    result = result * KT2KCAL # kcal/mol

                elif method == "cumu":
                    result_column_name = method + "_" + column_name
                    result = cumu(values)
                    result = result * KT2KCAL # kcal/mol
                else:
                    print("The method is not in method list!")
                    print("Abnormal Terminated")
                    sys.exit()

                results.append(result)


            df_result = pd.DataFrame()
            df_result[result_column_name] = results
            df_list.append(df_result)

    df_all = pd.concat(df_list, axis=1)
    df_all["FileNames"] = file_bases
    df_all.set_index("FileNames", inplace=True)
    df_all.to_csv("ani_d3_all_calculations.csv")


