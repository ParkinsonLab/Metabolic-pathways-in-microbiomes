#!/usr/bin/env python

import os
import os.path
import sys
import shutil

detect_ECs = sys.argv[1]
diamond_ECs = sys.argv[2]
priam_ECs = sys.argv[3]
con_dir = sys.argv[4]


try:
    os.mkdir(con_dir)
except:
    shutil.rmtree(con_dir)
    os.mkdir(con_dir)

with open(os.path.join(con_dir, "Consensus" + ".ECs_PD"), "w") as PB_out:
    with open(priam_ECs, "r") as priam_ECs_in:
        priam_preds = priam_ECs_in.readlines()
    with open(diamond_ECs, "r") as diamond_ECs_in:
        diamond_preds = diamond_ECs_in.readlines()
    PB_preds = []
    for priam_ec in priam_preds:
        if priam_ec in diamond_preds:
            PB_preds.append(priam_ec)
    PB_out.writelines(PB_preds)
with open(detect_ECs, "r") as detect_ECs_in:
    detect_preds = []
    for line in detect_ECs_in.readlines():
        line_as_list = line.split("\t")
        line_as_list = "\t".join(line_as_list[:2]) + "\n"
        detect_preds.append(line_as_list)
All_preds = set()
for pred in detect_preds:
    if len(pred.split("\t")[1].split(".")) == 4:
        All_preds.add(pred)
for pred in PB_preds:
    if len(pred.split("\t")[1].split(".")) == 4:
        All_preds.add(pred)
with open(os.path.join(con_dir, "Consensus" + ".ECs_All"), "w") as ec_out:
    ec_out.writelines(sorted(All_preds))