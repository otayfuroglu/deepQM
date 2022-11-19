#
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-csv_path", type=str, required=True, help="")
parser.add_argument("-start", type=int, required=False, default=0, help="")
parser.add_argument("-end", type=int, required=False, default=-1, help="")
parser.add_argument("-step", type=int, required=False, default=1, help="")

args = parser.parse_args()

csv_path = args.csv_path
start_point = args.start
end_point = args.end
step_size = args.step

df = pd.read_csv(csv_path)

list_target = []
for i, row in df.iloc[start_point:end_point].iterrows():
    if i % step_size == 0:
        list_target.append(row)

df = pd.concat(list_target, axis=1)
if end_point == -1:
    end_point = "end"
df.T.to_csv(f"from_{start_point}_to_{end_point}_by_{step_size}.csv", index=None)

# usege example
#python ../scripts/sampledWtihRange.py -csv_path 2A4F_fbf_energies.csv -start 0 -end 50 -step 5
