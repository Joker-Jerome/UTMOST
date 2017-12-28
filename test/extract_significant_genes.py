import numpy as np
import os 
import sys

outcome_file = sys.argv[1]
with open(outcome_file, "r") as fi:
	for line in fi:
		tmpline = line.strip().split()
		if tmpline[3] == "NA":
			continue
		if float(tmpline[3]) < 0.05 and float(tmpline[3]) > 0:
			print(tmpline[0] + "\n")

