#!/usr/local/bin/python2.7
# encoding: utf-8
'''
@author:     zhao_lab

@license:    license

@contact:    zhao_lab at Yale
'''


import os
import pandas as pd
import numpy as np

def match_list(a, b):
    """ find the indices of matching element 
    :param  a: list a 
            b: list b
    :return: The indices of elements in list a in list b 
    """
    return np.array([ b.index(x) if x in b else -1 for x in a ])

def run(args):
    tissue = args.output_file.split("/")[-1].split(".csv")[0]
    # filtering significant imputed genes
    filter_dir = os.path.dirname(os.path.abspath(args.model_db_path)) + "/glasso2_sig/"
    sig_gene_file = filter_dir + tissue + ".adj_expr.txt" 
    if os.path.exists(sig_gene_file):
        test = pd.read_csv(args.output_file)
        sig_gene = pd.read_csv(sig_gene_file, header = None)
        index = match_list(sig_gene.loc[:,0].tolist(), test["gene_name"].tolist())
        index_2 = np.setdiff1d(np.arange(0, len(test["gene_name"])), index[index > -1])
        test.loc[index_2, 'zscore'] = float('nan')
        test.loc[index_2, 'effect_size'] = float('nan')
        test.loc[index_2, 'pvalue'] = float('nan')
        test.to_csv(args.output_file, na_rep = "NA")
