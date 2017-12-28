import logging
import numpy
import pickle
import os
import gzip
import ntpath
import metax.WeightDBUtilities as WeightDBUtilities
import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.ThousandGenomesUtilities as ThousandGenomesUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.Formats as Formats
import sys
import sqlite3
import argparse

filename1 = sys.argv[1]
filename2 = sys.argv[2]
# file 1 
pickle_in = open(filename1,"rb")
tmp_logic_object = pickle.load(pickle_in)
pickle_in.close()

# file 2
pickle_in = open(filename2,"rb")
tmp_logic_object2 = pickle.load(pickle_in)
pickle_in.close()

print(len(tmp_logic_object.weights_by_gene.keys()))
print(len(tmp_logic_object2.weights_by_gene.keys()))
count = 0
gene_count = 0
diff_count = 0
for gene in tmp_logic_object.weights_by_gene.keys():
	count += len(tmp_logic_object.weights_by_gene[gene].keys())
	gene_count += 1
	print(gene)
	len1 = len(tmp_logic_object.weights_by_gene[gene].keys())
	len2 = len(tmp_logic_object2.weights_by_gene[gene].keys())
	if (len1!=len2):
		diff_count += 1
		print(len1)
print(count)	 
print(diff_count)
