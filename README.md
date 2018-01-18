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

* single_tissue_covariance.py

* single_tissue_association_test.py

* joint_tissue_covariance.py

* joint_GBJ_test.py

* test_tool

* metax module
                        
## Example and Usage

The following example assumes that you have **python 2.7**, **numpy**, **pandas**, **scipy**, **rpy2**, and **R** installed. 
All of these functions take different number of command line parameters. Run them with --help or -h option to see the options.


**1. Clone the UTMOST repository**
```bash
$ git clone https://github.com/Joker-Jerome/UTMOST
```

**2. Go to the software directory**
```bash
$ cd ./UTMOST
```

**3. Download sample data [sample data](dl.dropboxusercontent.com/s/pwk47cyiw60kkod/data.zip)**
```bash
# You can click on the link above or run the following
$ wget dl.dropboxusercontent.com/s/pwk47cyiw60kkod/data.zip?dl=0
```

**4. Unzip the data.zip file**
```bash
$ unzip data.zip
```
The data folder will include a **Model Database**, a **GWAS summary statistics**.

**5. Calculate the single tissue covariance**
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



**6. Run the single tissue association test**
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

* *--model_db_path* 

  Path to tissue transriptome model.
  
* *--covariance* 

  Path to file containing covariance information.
  
* *--gwas_folder* 

  Folder containing GWAS summary statistics data.
  
* *--gwas_file_pattern* 

  The file patten of gwas file.
  
* *--snp_column* 

  Argument with the name of the column containing the RSIDs.
  
* *--effect_allele_column* 

  Argument with the name of the column containing the effect allele.
  
* *--non_effect_allele_column* 

  Argument with the name of the column containing the non effect allele.
  
* *--beta_column* 

  The column containing -phenotype beta data for each SNP- in the input GWAS files.
  
* *--pvalue_column* 

  The column containing -PValue for each SNP- in the input GWAS files.
  
* *--output_file* 

  Path where results will be saved to.



**7. Calculate the joint tissue covariance**
```bash
python2 ./joint_tissue_covariance.py \
--weight_db /ysm-gpfs/home/zy92/project/UTMOST/database/weigth_db_v2/ \
--input_folder /ysm-gpfs/home/zy92/project/UTMOST/database/dosage/ \
--covariance_output /ysm-gpfs/home/zy92/project/UTMOST/intermediate/

```

The example command parameters:

* *--verbosity* 
  
  Log verbosity level. 1 means everything will be logged. 10 means high level messages will be logged.
  
* *--weight_db*

  Name of weight db in data folder.
  
* *--input_folder*

  Name of folder containing dosage data.
  
* *--covariance_output* 

  Path where covariance results will be saved to.
  
* *--min_maf_filter*

  Filter snps according to this maf.
  
* *--max_maf_filter*

  Filter snps according to this maf.
  

8. Joint GBJ test
```bash
$ python2 joint_GBJ_test.py \
--weight_db /ysm-gpfs/home/zy92/project/UTMOST/database/weigth_db_v2/ \
--output_dir /ysm-gpfs/home/zy92/project/UTMOST/outcome/ \
--cov_dir /ysm-gpfs/home/zy92/project/UTMOST/intermediate/ \
--input_folder /ysm-gpfs/home/zy92/project/UTMOST/database/mask/AD/ \
--gene_info /ysm-gpfs/home/zy92/project/UTMOST/intermediate/gene_info.txt \
--output_name test
```

The example command parameters:

* *--verbosity* 

  Log verbosity level. 1 means everything will be logged. 10 means high level messages will be logged.
  
* *--weight_db*

  Name of weight db in data folder.
  
* *--input_folder*

  Name of folder containing summary data.
  
* *--cov_dir* 

  Path where covariance results are.
  
* *--output_dir* 

  Path where results will be saved to.
  
* *--gene_info*

  Name of file containing the gene list.
  
* *--start_gene_index*

  Index of the starting gene.
  
* *--end_gene_index*

  Index of the ending gene
  
## Acknowledgement

## Reference
