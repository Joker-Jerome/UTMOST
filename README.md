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

All of these functions take different number of command line parameters. Run them with --help or -h option to see the options.

### Joint Covariance Calculation

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
                        
## Example and Usage
The following example assumes that you have **python 2.7**, **numpy**, **rpy2**, **scipy** and **R** installed.

1. Clone the UTMOST repository 
```bash
$ git clone https://github.com/Joker-Jerome/UTMOST
```

2. Go to the software directory
```bash
$ cd ./UTMOST
```

3. Download sample data [sample data](dl.dropboxusercontent.com/s/pwk47cyiw60kkod/data.zip):
```bash
# You can click on the link above or run the following
$ wget dl.dropboxusercontent.com/s/pwk47cyiw60kkod/data.zip?dl=0
```

4. Unzip the data.zip file
```bash
$ unzip data.zip
```
The data folder will include a **Model Database**, a **GWAS summary statistics**.

5. Calculate the single tissue covariance
```bash
python2 ./single_tissue_covariance.py \
--weight_db data/DGN-WB_0.5.db \
--input_folder data/GWAS/ \
--covariance_output data/covariance.txt.gz
```
The example command parameters:

* *--model_db_path* Path to tissue transriptome model.
* *--input_folder* Folder containing GWAS summary statistics data.
* *--covariance_output* Path where covariance will be saved to.

6. Run the single tissue association test
```bash
python2 ./single_tissue_association_test.py \
--model_db_path data/DGN-WB_0.5.db \
--covariance data/covariance.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/single_tissue_test_results.csv
```
The example command parameters:

* *--model_db_path* Path to tissue transriptome model.
* *--covariance* Path to file containing covariance information.
* *--gwas_folder* Folder containing GWAS summary statistics data.
* *--gwas_file_pattern* The file patten of gwas file.
* *--snp_column* Argument with the name of the column containing the RSIDs.
* *--effect_allele_column* Argument with the name of the column containing the effect allele.
* *--non_effect_allele_column* Argument with the name of the column containing the non effect allele.
* *--beta_column* The column containing -phenotype beta data for each SNP- in the input GWAS files.
* *--pvalue_column* The column containing -PValue for each SNP- in the input GWAS files.
* *--output_file* Path where results will be saved to.
  
## Acknowledgement

## Reference
