import os

os.chdir("/home/chase/DissertationProjects/nextstrain/WNV_Nextstrain_Analysis/")

with open("nextstrain/data/metadata.tsv", "r") as meta_file:

    while line := meta_file.readline():
        print(line)

