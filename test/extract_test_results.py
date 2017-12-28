import numpy as np
import os 
import sys

outcome_file = sys.argv[1]
test_dict = {}
with open(outcome_file, "r") as fi:
	for line in fi:
		tmpline = line.strip().split()
		gene = tmpline[0]
		if gene not in test_dict.keys():
			test_dict[gene] = [tmpline[2],tmpline[3]]
		if tmpline[3] == "NA":
			continue
		if float(tmpline[3]) < 0.05 and float(tmpline[3]) > 0:
			print(tmpline[0] + "\n")

output_file = sys.argv[2]
gene_list = sorted(list(test_dict.keys()))
with open(output_file,"w") as fo:
	for gene in gene_list:
		tmpline = "\t".join([gene,test_dict[gene][0],test_dict[gene][1]])
		fo.write(tmpline + "\n")


