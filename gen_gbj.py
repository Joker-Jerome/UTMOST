'''
Created on Dec 30, 2017

@author: jerome
'''

import os
import sys

disease = sys.argv[1]
split = int(sys.argv[2])

num = 15605
step = (num - (num % split))/split

for i in range(split):
    tmp_filename = disease + "_" + "split" + str(i+1) + ".src"
    with open(tmp_filename,"w") as fo:
        tmp_start = int(i * step + 1)
        if i == split - 1:  
            tmp_end = num
        else:    
            tmp_end = int((i + 1) * step) 
        fo.write("#!/bin/bash\n")
        fo.write("#SBATCH -p general\n")
        fo.write("#SBATCH -J " + disease + "_gbj\n")
        fo.write("#SBATCH -n 1 --cpus-per-task=6\n")
        fo.write("#SBATCH --mem-per-cpu=12G\n")
        fo.write("date\n")
        fo.write("python joint_GBJ.py --weight_db ~/scratch60/database_tissue/ --input_folder /ysm-gpfs/pi/zhao/ml2376/association_v3/" + disease\
                 + "/single_mask/" + " --start_gene_index " + str(tmp_start) + " --end_gene_index " + str(tmp_end) + " --output_name "  + disease)
        fo.write("date\n")
