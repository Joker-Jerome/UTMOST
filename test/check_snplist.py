import numpy as np
import os 
import pandas as pd
import sys

gene_list = pd.read_table("/ysm-gpfs/home/zy92/project/metaxcan/MetaXcan/software/gene_info.txt").iloc[:,0].values
dir_a = sys.argv[1]
dir_b = sys.argv[2]
diff_list = []
more_list = []
less_list = []
exist_list = []
print(gene_list[:5])
count = 0
for gene in gene_list:
	count += 1
	if count % 1000 == 0:
		print(str(count) + " genes have been processed.\n")
	filename_a = dir_a + "/" + gene + ".snplist"
	filename_b = dir_b + "/" + gene + ".snplist"
	if os.path.exists(filename_a) and os.path.exists(filename_b):
		list_a = pd.read_table(filename_a, header = None).iloc[:,0].values
		list_b = pd.read_table(filename_b, header = None).iloc[:,0].values
		if (set(list_a.tolist()) != set(list_b.tolist())):
			diff_list.append(gene)
		if (len(list_a) > len(list_b)):
			more_list.append(gene)
			print(gene)
			print("list_a:")
			print(len(list_a))
			print("list_b:")
			print(len(list_b))
		elif (len(list_a) < len(list_b)):
			less_list.append(gene)
			print(gene)
			print("list_a:")
			print(len(list_a))
			print("list_b:")
			print(len(list_b))
	else:
		exist_list.append(gene)
print("diff_list")	
print(len(diff_list))	
print("more_list")
print(len(more_list))
print("less_list")
print(len(less_list))
print("exist_list")
print(len(exist_list))

