#!/usr/bin/env python

import sys

Input_File = sys.argv[1]
Output_File = sys.argv[2]

with open(Input_File, "r") as ECs:
    with open(Output_File, "w") as processedECs:
        for line in ECs.readlines():
            columns = line.split("\t")
            if len(columns[1].split(";")) > 1:
                for EC in range(len(columns[1].split(";"))):
                    processedECs.write("\t".join([columns[0], columns[1].split(";")[EC].strip()]))
                    processedECs.write("\n")
            else:
                processedECs.write(line)