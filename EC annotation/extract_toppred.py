#!/usr/bin/env python

import sys

Input_File = sys.argv[1]
Output_File = sys.argv[2]

with open(Input_File, "r") as topred:
    with open(Output_File, "w") as cutoff:
        for line in topred:
            columns = line.split("\t")
            if columns[2] == "probability":
                continue
            if float(columns[2]) >= 0.2 and int(columns[3]) > 5:
                cutoff.write(line)