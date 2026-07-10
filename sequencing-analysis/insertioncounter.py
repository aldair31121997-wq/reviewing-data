#!/usr/bin/env python
# coding: utf-8

## load modules
import pandas as pd
import argparse
from collections import Counter

# load dataframe
parser = argparse.ArgumentParser(description="Process compiled mutpos table")
parser.add_argument("input_table", help="Input TSV file, e.g. compilatedclean.tab")

args = parser.parse_args()

df = pd.read_table(args.input_table, sep="\t")

# change the data type of the position.
df = df.astype({"Pos": "int"})

# mask

# get the replica column
df["rep"] = df["sample"].str.strip().str[-1]

# remove replica information from sample name
df["sample"] = (
    df["sample"]
    .str.replace("-R1", "")
    .str.replace("-R2", "")
    .str.replace("-R3", "")
    .str.replace("-R4", "")
    .str.replace("-R5", "")
)

df["sample"].unique()

df["cell"] = df["sample"].str.split("-", n=1).str[1]

df["sample"] = df["sample"].str.split("-").str[0]

# create freqmuts column for further analysis
df["FreqMuts"] = df["Muts"] / df["Depths"]
df["FreqInser"] = df["inscount"] / df["Depths"]
df["FreqDel"] = df["delcount"] / df["Depths"]

# getme the indentity of all insertion observed
mask_inser = (df["Pos"] >= 1000) & (df["Pos"] <= 1050)

dfinser = (
    df.loc[mask_inser]
    .groupby(
        [
            "sample",
            "Pos",
            "inscount",
            "FreqInser",
            "Template",
            "Depths",
            "rep",
            "cell",
        ]
    )["ins_seq_list"]
    .apply(list)
)

dfinser = pd.DataFrame(dfinser)
dfinser.reset_index(inplace=True)

completeinser2 = pd.DataFrame(
    columns=[
        "sample",
        "Pos",
        "maxrate",
        "insertion",
        "FreqInser",
        "Template",
        "Depths",
        "rep",
        "cell",
    ]
)

counter = 0

# eliminate NAS
for r in range(len(dfinser["sample"])):
    dfinser.at[r, "ins_seq_list"] = [
        x for x in dfinser["ins_seq_list"][r] if str(x) != "nan"
    ]

    # this long line takes the elements on the insertion list and evaluates their frequence.
    # the above for allows to iterate over the whole dataframe
    d = Counter(
        "".join(dfinser["ins_seq_list"][r])
        .replace("'", "")
        .replace("[", "")
        .replace("]", "")
        .replace(" ", "")
        .replace(",", "")
        .replace("nan", "")
        .split("+")
    )

    j = d.most_common(1)[0][0]
    test = d.most_common(3)

    insertions = pd.DataFrame.from_dict(d, orient="index").reset_index()
    insertions["Pos"] = dfinser["Pos"][r]
    insertions["sample"] = dfinser["sample"][r]

    for p in range(len(test)):
        completeinser2.loc[counter] = (
            [dfinser["sample"][r]]
            + [dfinser["Pos"][r]]
            + [test[p][1]]
            + [test[p][0]]
            + [dfinser["FreqInser"][r]]
            + [dfinser["Template"][r]]
            + [dfinser["Depths"][r]]
            + [dfinser["rep"][r]]
            + [dfinser["cell"][r]]
        )

        counter = counter + 1

completeinser2["FreqInser"] = completeinser2["maxrate"] / completeinser2["Depths"]

# export to excel format for further analysis
completeinser2.to_excel("insertions.xlsx")
