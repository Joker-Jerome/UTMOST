# UTMOST

UTMOST is a tool to perform integrative genetic studies. 

## Prerequisites

The software is developed and tested in Linux and Mac OS environments. 

To run UTMOST, you need Python 2.7, numpy (>=1.11.1), scipy (>=0.18.1), pandas (>=0.18.1), rpy2 (==2.8.6).

R is needed for GBJ testing.

## Project Layout

## Input data

## Output Format

## Usage

### Joint GBJ Test
joint_GBJ.py

usage: joint_GBJ.py 
[-h] 
[--verbosity VERBOSITY] 
[--weight_db WEIGHT_DB]
[--output_dir OUTPUT_DIR] 
[--cov_dir COV_DIR]
[--input_folder INPUT_FOLDER] 
[--gene_info GENE_INFO]                    
[--start_gene_index START_GENE_INDEX]
[--end_gene_index END_GENE_INDEX]

optional arguments:

  -h, --help            show this help message and exit
  --verbosity VERBOSITY
                        Log verbosity level. 1 is everything being logged. 10
                        is only high level messages, above 10 will hardly log
                        anything
  --weight_db WEIGHT_DB
                        name of weight db in data folder
  --output_dir OUTPUT_DIR
                        the output directory
  --cov_dir COV_DIR     the covariance directory
  --input_folder INPUT_FOLDER
                        name of folder containing summary data
  --gene_info GENE_INFO
                        name of folder containing summary data
  --start_gene_index START_GENE_INDEX
                        index of the starting gene
  --end_gene_index END_GENE_INDEX
                        index of the ending gene

## Acknowledgement

## Reference
