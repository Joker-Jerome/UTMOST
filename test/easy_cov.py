'''
Created on Dec 14, 2017

@author: jerome
'''
import sys
import os
import numpy as np
import pandas as pd
import gzip

rsid1 = sys.argv[1]
rsid2 = sys.argv[2]
seq1 = []
seq2 = []

dosage_path = "/ysm-gpfs/pi/zhao/ml2376/1000g_phase3/dosage/" 

fi = []
for file in sorted(os.listdir(dosage_path)):
    if file.endswith("dosage.gz"):
        fi.append(file)
print(fi)
print(rsid1)
print(rsid2)
os.chdir(dosage_path)

for dosage_file in fi:
    with gzip.open(dosage_file,"rb") as fi:
        for line in fi:
            tmpline = line.strip().split()
            if tmpline[1].decode("utf-8") == rsid1:
                seq1 = list(map(lambda x: float(x.decode("utf-8")),tmpline[6:]))
                print(seq1)
            if tmpline[1].decode("utf-8") == rsid2:
                seq2 = list(map(lambda x: float(x.decode("utf-8")),tmpline[6:]))
                print(seq2)
            if (len(seq1) != 0) and (len(seq2) != 0):
                print(np.cov(seq1, seq2))   
                break

 

