#
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-csv_path", type=str, required=True, help="")

args = parser.parse_args()
csv_path = args.csv_path

def getBlocks(csv_path: str) -> dict:
    flag = False
    blocks_dict = {}
    with open(csv_path) as lines:
        for i, line in enumerate(lines):
            if ":" in line:
                fl_name = ("%s.csv"%line
                           .replace("\n", "")
                           .replace(":", "")
                           .replace(" ", "_")
                           .lower())
                blocks_dict[fl_name] = []
                #  fl_block = open(fl_name, "w")
                flag = True
            if flag:
                blocks_dict[fl_name].append(line)
    return blocks_dict


def getSplitMethodBlock(blocks_dict):
    for fl_name, lines in blocks_dict.items():
        block_dict = {}
        block = []
        key = 0
        #  with open(csv_path) as lines:
        for i, line in enumerate(lines):
            block += [line]
            if "Energy Terms" in line:
                block_dict[str(key)] = block
                block = []
                key += 1
        # 0. item is not neccessary
        del block_dict["0"]

        fl_block = open(fl_name, "w")
        for j in range(len(block_dict["1"])):
            for i in range(len(block_dict.keys())):
                line = block_dict[str(i+1)][j]
                if "Energy Terms" in line:
                    continue
                if i == 2:
                    fl_block.write(line)
                else:
                    fl_block.write(line.replace("\n", "")+",")
        fl_block.close()

# usage example
#python ../scripts/sampledWtihRange.py -csv_path example.csv 
if __name__ == "__main__":
    blocks_dict = getBlocks(csv_path)
    getSplitMethodBlock(blocks_dict)
