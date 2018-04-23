# UTMOST

UTMOST (Unified Test for MOlecular SignaTures) is a principled method to perform cross-tissue expression imputation and gene-level association analysis. The preprint could be found at [A statistical framework for cross-tissue transcriptome-wide association analysis](https://www.biorxiv.org/content/early/2018/03/21/286013).

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

* joint_covariance.py

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
**3. Download reference data (~12GB for zipped file, ~50GB after unzipping)**
```bash
# You can click on the link above or run the following
$ wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies  /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1u8CRwb6rZ-gSPl89qm3tKpJArUT8XrEe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1u8CRwb6rZ-gSPl89qm3tKpJArUT8XrEe" -O sample_data.zip && rm -rf /tmp/cookies.txt
$ unzip sample_data.zip
```
The data folder will include the following files/folders:
```
weight_db_v2
```
**4. Quick start**
**The following part is a demonstration of applying UTMOST with imputation models jointly trained in 44 tissues with GTEx data. The sample_data.zip contains pre-calculated imputation models, covariance matrices for single tissue and joint tissue GBJ test. Pipeline for generating your own covariance matrices and incorporating other eQTL/sQTL data (e.g. from [STARNET](https://github.com/Wainberg/Vulnerabilities_of_TWAS) and [BLUEPRINT]((ftp://ftp.ebi.ac.uk/pub/databases/blueprint/blueprint Epivar/qtl_as/), i.e. the analysis pipeline used in manuscript) could be found in Section 5.**

**4.1. Run single tissue association test for 44 tissues**
```bash
TISSUE_GTEx=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis, Thyroid Uterus Vagina Whole_Blood)
for tissue in ${TISSUE_GTEx[@]}
do
python2 ./single_tissue_association_test.py \
--model_db_path sample_data/weight_db_v2/${tissue}.db \
--covariance sample_data/covariance/${tissue}.txt.gz \
--gwas_folder sample_data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/${tissue}.csv
done
```
The example command parameters:

* *--model_db_path* 

  Path to gene expression imputation model (estimated weights/effect sizes of cis-eQTLs).
  
* *--covariance* 

  Path to file containing covariance information (used to estimate the variance of gene-level effect size estimator, see **Gene-level association test** in **Methods** section of manuscript).
  
* *--gwas_folder* 

  Folder containing GWAS summary statistics data.
  
* *--gwas_file_pattern* 

  The file patten of gwas file (file name of summary statistics if not segmented by chromosomes).
  
* *--snp_column* 

  Argument with the name of the column containing the RSIDs.
  
* *--effect_allele_column* 

  Argument with the name of the column containing the effect allele.
  
* *--non_effect_allele_column* 

  Argument with the name of the column containing the non-effect allele.
  
* *--beta_column* 

  The column containing -effect size estimator for each SNP- in the input GWAS files.
  
* *--pvalue_column* 

  The column containing -PValue for each SNP- in the input GWAS files.
  
* *--output_file* 

  Path where results will be saved to.

**4.2. Combine gene-trait associations in 44 tissues by joint GBJ test**
```bash
mkdir sample_data/results_joint ## save association results for cross-tissue joint test
$ python2 joint_GBJ_test.py \
--weight_db sample_data/weight_db_v2/ \
--output_dir sample_data/results_joint/ \
--cov_dir sample_data/covariance_joint/ \
--input_folder sample_data/results/ \
--gene_info intermediate/gene_info.txt \
--output_name test
```

The example command parameters:

* *--verbosity* 

  Log verbosity level. 1 means everything will be logged. 10 means high level messages will be logged.
  
* *--weight_db*

  Name of weight db in data folder (imputation models).
  
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

Output format:

  | Gene        | Test score     | P value  |
  | ------------- |:-------------:| -----:|
  |      Gene A    | test score A               |  P value A        |
  |      Gene B    | test score B               |  P value B        |

**5. Calculate the single tissue covariance**
```bash
python2 ./single_tissue_covariance.py \
--weight_db sample_data/DGN-WB_0.5.db \
--input_folder sample_data/dosage/ \
--covariance_output sample_data/covariance.txt.gz
```
The example command parameters:

* *--weight_db* 

  Path to tissue transriptome model.
  
* *--input_folder* 

  Folder containing GWAS summary statistics data.
  
* *--covariance_output* 

  Path where covariance will be saved to.



**6. Run the single tissue association test**
```bash
python2 ./single_tissue_association_test.py \
--model_db_path sample_data/DGN-WB_0.5.db \
--covariance sample_data/covariance.txt.gz \
--gwas_folder sample_data/GWAS \
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
python2 ./joint_covariance.py \
--weight_db /YOUR_UTMOST_DIR/sample_data/weight_db_v2/ \
--input_folder /YOUR_UTMOST_DIR/YOUR_DOSAGE_DIR/ \
--covariance_output /YOUR_UTMOST_DIR/intermediate/tmp/

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
--weight_db /YOUR_UTMOST_DIR/sample_data/weight_db_v2/ \
--output_dir /YOUR_UTMOST_DIR/results/ \
--cov_dir /YOUR_UTMOST_DIR/sample_data/covariance_joint/ \
--input_folder /YOUR_UTMOST_DIR/sample_data/mask/ \
--gene_info /YOUR_UTMOST_DIR/intermediate/gene_info.txt \
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

Output format:

  | Gene        | Test score     | P value  |
  | ------------- |:-------------:| -----:|
  |      Gene A    | test score A               |  P value A        |
  |      Gene B    | test score B               |  P value B        |

  
## Acknowledgement
Part of the code is modified from MetaXcan https://github.com/hakyimlab/MetaXcan. We thank them for sharing the code.

## Reference
**Hu et al. (2018). A statistical framework for cross-tissue transcriptome-wide association analysis. bioRxiv, 286013.**
[Link](https://www.biorxiv.org/content/early/2018/03/21/286013)
