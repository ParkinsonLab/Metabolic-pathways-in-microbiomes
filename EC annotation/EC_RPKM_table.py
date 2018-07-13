#!/usr/bin/env python

import sys

Gene2EC = sys.argv[1]
Gene_table = sys.argv[2]
EC_table = sys.argv[3]

EC2genes = {}
with open(Gene2EC, "r") as infile:
    for line in infile:
        columns = line.split("\t")
        gene = columns[0]
        EC = columns[1].strip("\n")
        if EC in EC2genes:
            EC2genes[EC].append(gene)
        else:
            EC2genes[EC] = [gene]

Gene_dict = {}
Header = True
Header_text = ""
with open(Gene_table, "r") as infile:
    for line in Gene_table:
        if Header:
            Header_text = line.strip("\n")
            Header = False
        else:
            columns = line.split(",")
            gene = columns[0]
            RPKMs = []
            for RPKM in columns[1:]:
                RPKMs.append(RPKM.strip("\n"))
            Gene_dict[gene] = RPKMs

EC_dict = {}
for EC in EC2genes:
    EC_dict[EC] = []
    for gene in EC2genes[EC]:
        if EC_dict[EC] == []:
            EC_dict[EC] = Gene_dict[gene]
        else:
            RPKMs = []
            for n in range(len(EC_dict[EC])):
                RPKMs.append(str(float(EC_dict[EC][n]) + float(Gene_dict[gene][n])))
            EC_dict[EC] = RPKMs

with open(EC_table, "w") as outfile:
    outfile.write(Header_text + "\n")
    for EC in EC_dict:
        outfile.write(EC + "," + ",".join(EC_dict[EC]) + "\n")