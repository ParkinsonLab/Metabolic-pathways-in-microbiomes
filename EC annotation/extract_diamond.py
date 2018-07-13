#!/usr/bin/env python

import sys

Input_File = sys.argv[1]
Output_File = sys.argv[2]
SWISS_PROT_MAP = sys.argv[4]

mapping_dict = {}

with open(SWISS_PROT_MAP, "r") as mapping:
    for line in mapping.readlines():
        columns = line.split("\t")
        mapping_dict[columns[0]] = set(columns[2:])

with open(Input_File, "r") as blastout:
    with open(Output_File, "w") as ecout:
        for line in blastout.readlines():
            columns = line.strip().split("\t")
            for EC in mapping_dict:
                if columns[1] in mapping_dict[EC]:
                    ecout.write("\t".join([columns[0], EC + "\n"]))