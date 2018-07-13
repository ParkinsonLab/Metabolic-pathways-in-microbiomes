#!/usr/bin/env python

import os
import os.path
import sys

splitlen = int(sys.argv[1])
input_file = sys.argv[2]
output_folder = sys.argv[3]

def write_fasta(fasta_dict, count, mode):
    output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0] + "_split_" + str(count).zfill(3) + ".fasta")
    with open(output_file, mode) as outfile:
        for seq_id in fasta_dict:
            outfile.write(">" + seq_id + "\n")
            outfile.write(fasta_dict[seq_id][0] + "\n")

os.chdir(os.path.dirname(input_file))

prot_count = 0
with open(input_file, "r") as infile:
    for line in infile:
        if line.startswith(">"):
            prot_count += 1
File_count = prot_count / splitlen + 1
for n in range(File_count):
    Start = n*splitlen
    m = n + 1
    Stop = m*splitlen
    Temp_sequences = {}
    with open(input_file, "r") as infile:
        seq_id = ""
        seq_count = 0
        for line in infile:
            if line.startswith(">"):
                if seq_count < Start:
                    seq_count += 1
                elif seq_count >= Start and seq_count < Stop:
                    seq_id = line[1:].split(" ")[0].strip("\n")
                    Temp_sequences[seq_id] = [""]
                    seq_count += 1
                elif seq_count == Stop:
                    break
            elif seq_id != "":
                Temp_sequences[seq_id][0] += line.strip("\n")
    write_fasta(Temp_sequences, m, "w")