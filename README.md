# UTMOST

UTMOST is a tool to perform integrative genetic studies. 

## Prerequisites

The software is developed and tested in Linux and Mac OS environments. 

* Python 2.7 

* numpy (>=1.11.1)

* scipy (>=0.18.1)

* pandas (>=0.18.1)

* rpy2 (==2.8.6)

* R is needed for GBJ testing.

## Project Layout

* Covariance.py

* joint_GBJ.py

* test_tool

* metax module


## Usage

### Covariance Calculation

* usage: Covariance.py [-h] 

[--verbosity VERBOSITY]
                              
[--weight_db WEIGHT_DB]

[--input_folder INPUT_FOLDER]

[--covariance_output COVARIANCE_OUTPUT]

[--input_format INPUT_FORMAT]

[--min_maf_filter MIN_MAF_FILTER]

[--max_maf_filter MAX_MAF_FILTER]

[--max_snps_in_gene MAX_SNPS_IN_GENE]

Build covariances from dosage data and weights database.

optional arguments:
  -h, --help            Show this help message and exit
  
  --verbosity VERBOSITY
  
                        Log verbosity level. 1 is everything being logged. 10
                        is only high level messages, above 10 will hardly log
                        anything
  --weight_db WEIGHT_DB
  
                        Name of weight db in data folder
                        
  --input_folder INPUT_FOLDER
  
                        Name of folder containing dosage data
                        
  --covariance_output COVARIANCE_OUTPUT
  
                        Name of file to dump covariance results in. Defaults
                        to 'intermediate/cov/' + file name prefix from '--
                        weight_db' argument
                        
  --input_format INPUT_FORMAT
  
                        Input dosage files format. Valid options are: IMPUTE,
                        PrediXcan
                        
  --min_maf_filter MIN_MAF_FILTER
  
                        Filter snps according to this maf
                        
  --max_maf_filter MAX_MAF_FILTER
  
                        Filter snps according to this maf
                        
  --max_snps_in_gene MAX_SNPS_IN_GENE
  
                        Ignore any gene that has snps above this value
### Joint GBJ Test

* usage: joint_GBJ.py 

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
  
                        Name of weight db in data folder
                        
  --output_dir OUTPUT_DIR
  
                        Output directory
                        
  --cov_dir COV_DIR     
  
                        Covariance directory
  
  --input_folder INPUT_FOLDER
  
                        Name of folder containing summary data
                        
  --gene_info GENE_INFO
  
                        Name of folder containing summary data
                        
  --start_gene_index START_GENE_INDEX
  
                        Index of the starting gene
                        
  --end_gene_index END_GENE_INDEX
  
                        Index of the ending gene
                        

## Acknowledgement

## Reference
