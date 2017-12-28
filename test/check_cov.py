'''
Created on Dec 13, 2017

@author: jerome
'''
import sys
import os
import numpy as np
import pandas as pd

res_path = sys.argv[1]
gene = sys.argv[2]

snp_file = res_path + gene + ".snplist"
cov_file = res_path + gene + ".cov"

snp_rsid = pd.read_table(snp_file, header = None)
snp_list = list(snp_rsid.loc[:,0])
print("Number of Snps: " + str(len(snp_list)))
sorted_snp_list = sorted(snp_list)

n = len(snp_list)

cov_matrix = np.loadtxt(cov_file)

for i in range(n):
    for j in range(n):
        cov_i = snp_list.index(sorted_snp_list[i])
        cov_j = snp_list.index(sorted_snp_list[j])
        print(sorted_snp_list[i] + "\t" + sorted_snp_list[j] + "\t" + str(cov_matrix[cov_i,cov_j]))


