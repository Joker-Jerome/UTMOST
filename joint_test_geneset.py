#!/usr/local/bin/python2.7
# encoding: utf-8
'''
joint_GBJ_geneset

@contact:    zhao_lab at Yale

'''

import pickle
import sys
import os
import gzip
import ntpath
import metax
import logging
import readline
import sqlite3
import re
import subprocess
import numpy as np
import pandas as pd
import pandas 
import copy
import metax.WeightDBUtilities as WeightDBUtilities
import metax.PrediXcanFormatUtilities as PrediXcanFormatUtilities
import metax.ThousandGenomesUtilities as ThousandGenomesUtilities
import metax.Logging as Logging
import metax.Utilities as Utilities
import metax.Formats as Formats
import metax.DataSetSNP as DataSetSNP
import scipy.stats as st


from timeit import default_timer as timer
from metax import Constants
from metax.gwas import GWAS
from metax.gwas import Utilities as GWASUtilities
from metax import PredictionModel
from metax import Utilities
from metax import Logging
from metax import Exceptions


def pathLeaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def mergeTwoDicts(x, y):
    z = x.copy()  
    z.update(y)   
    return z

# function of creating connection
def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except sqlite3.Error as e:
        print(e)
        sys.exit(1)
    return None

# function of matching two lists
def match_list(a, b):
    """ find the indices of matching element 
    :param  a: list a 
            b: list b
    :return: The indices of elements in list a in list b 
    """
    return np.array([ b.index(x) if x in b else -1 for x in a ])


def align_data_to_alleles(data, base, left_on, right_on):
    EA, NEA = Constants.EFFECT_ALLELE, Constants.NON_EFFECT_ALLELE
    EA_BASE, NEA_BASE = EA+"_BASE", NEA+"_BASE"
    merged = pandas.merge(data, base, left_on=left_on, right_on=right_on, suffixes=("", "_BASE"))

    alleles_1 = pandas.Series([set(e) for e in zip(merged[EA], merged[NEA])])
    alleles_2 = pandas.Series([set(e) for e in zip(merged[EA_BASE], merged[NEA_BASE])])
    eq = alleles_1 == alleles_2
    merged = merged[eq]

    flipped = merged[EA] != merged[EA_BASE]
    Z = Constants.ZSCORE
    if Z in merged:
        merged.loc[flipped, Z] = - merged.loc[flipped, Z]
    B = Constants.BETA
    if B in merged:
        merged.loc[flipped, B] = - merged.loc[flipped, B]

    merged.loc[flipped, EA] = merged.loc[flipped, EA_BASE]
    merged.loc[flipped, NEA] = merged.loc[flipped, NEA_BASE]

    return merged

def build_betas(args, model, gwas_format, name):
    load_from = os.path.join(args.gwas_folder, name)
    if model or args.skip_until_header:
        snps = model.snps() if model else None
        snp_column_name = args.snp_column if model else None
        load_from = GWASUtilities.gwas_filtered_source(load_from, snps=snps, snp_column_name=snp_column_name, skip_until_header=args.skip_until_header, separator=args.separator)
    sep = '\s+' if args.separator is None else args.separator
    b = GWAS.load_gwas(load_from, gwas_format, sep=sep, input_pvalue_fix=args.input_pvalue_fix)

    if model is not None:
        PF = PredictionModel.WDBQF
        base = model.weights[[PF.K_RSID, PF.K_EFFECT_ALLELE, PF.K_NON_EFFECT_ALLELE]].drop_duplicates()
        b = align_data_to_alleles(b, base, Constants.SNP, PF.K_RSID)

    b = b.fillna("NA")
    keep = [GWAS.SNP, GWAS.EFFECT_ALLELE, GWAS.NON_EFFECT_ALLELE, GWAS.ZSCORE]
    if GWAS.BETA in b: keep.append(GWAS.BETA)
    b = b[keep]
    return b

def validate(args):
    if not args.gwas_folder: raise Exceptions.InvalidArguments("You need to provide an input folder containing GWAS files")

def readGWAS(args):
    start = timer()
    validate(args)
    regexp = re.compile(args.gwas_file_pattern) if args.gwas_file_pattern else  None
    names = Utilities.contentsWithRegexpFromFolder(args.gwas_folder, regexp)
    names.sort() #cosmetic, because different filesystems/OS yield folders in different order

    if len(names) == 0:
        msg = "No GWAS files found on %s with pattern %s" % (args.gwas_folder, args.gwas_file_pattern,)
        raise Exceptions.ReportableException(msg)
    
    print "INFO: Reading GWAS data"
    gwas_format = GWASUtilities.gwas_format_from_args(args)
    GWAS.validate_format_basic(gwas_format)
    GWAS.validate_format_for_strict(gwas_format)
    #model = PredictionModel.load_model(args.model_db_path) if args.model_db_path else None
    model = None
    # dataframe
    r = pandas.DataFrame()
    for name in names:
        b = build_betas(args, model, gwas_format, name)
        r = pandas.concat([r,b])
    end = timer()
    logging.info("Successfully parsed input gwas in %s seconds"%(str(end-start)))
    print("Successfully parsed input gwas in %s seconds"%(str(end-start)))
    return r

def buildGWAS(dataframe, snplist, output_path, flip_tag_vec):
    tmp_df = copy.deepcopy(dataframe[dataframe['snp'].isin(snplist)])
    snplist = tmp_df['snp'].values
    snp2tag = dict(zip(snplist, flip_tag_vec))
    flip_tag_vec = map(lambda x:snp2tag[x], tmp_df['snp'].values)
    tmp_df['zscore'] = tmp_df['zscore'].values * np.array(flip_tag_vec)
    tmp_df['beta'] = tmp_df['beta'].values * np.array(flip_tag_vec)
    tmp_df.to_csv(output_path, sep='\t',index=None)
        
    
class PDTF:
    """Format of PrediXcan dosage"""
    CHR = 0
    RSID = 1
    POSITION = 2
    ALLELE_0 = 3
    ALLELE_1 = 4
    COLUMN_5 = 5
    FIRST_DATA_COLUMN = 6
    
    
class PrediXcanFormatDosageLoader(object):
    #weight_db_logic is used only for discarding absent snps.
    def __init__(self, path, weight_db_logic):
        self.path = path
        self.weight_db_logic = weight_db_logic

    def load(self):
        #print "INFO: Loading dosage files"
        #logging.info("Loading %s dosage", self.path)
        class PrediXcanCollector(object):
            def __init__(self, snps=[], snps_by_rsid={}, weight_db_logic=None):
                self.snps = snps
                self.snps_by_rsid = snps_by_rsid
                self.weight_db_logic = weight_db_logic

            def __call__(self, i, components):
                rsid = components[PDTF.RSID]
                if self.weight_db_logic and not rsid in self.weight_db_logic.genes_for_an_rsid:
                    logging.log(5, "rsid %s not in weight db, skip it", rsid)
                    return

                position = components[PDTF.POSITION]

                ref_allele = components[PDTF.ALLELE_0]
                if not ref_allele in Utilities.VALID_ALLELES:
                    logging.log(9, "wrong ref allele, rsid %s is not an SNP", rsid)
                    return
                eff_allele = components[PDTF.ALLELE_1]
                if not eff_allele in Utilities.VALID_ALLELES:
                    logging.log(9, "wrong eff allele, rsid %s is not an SNP", rsid)
                    return
                dosages = map(float,components[PDTF.FIRST_DATA_COLUMN:]) #dosages may be inputed
                #Should we flip based on weight_db at this point?

                snp = DataSetSNP.DataSetSNP(name=rsid, index=i, data=dosages, position=int(position), ref_allele=ref_allele, eff_allele=eff_allele)
                if snp.name in self.snps_by_rsid:
                    old = self.snps_by_rsid[snp.name]
                    logging.info("Duplicated rsid: (%s,%s) %s", old.name, old.position, " ".join(components))
                self.snps.append(snp)
                self.snps_by_rsid[snp.name] = snp
        loader = Utilities.CSVFileIterator(self.path, compressed=True)
        collector = PrediXcanCollector(weight_db_logic=self.weight_db_logic)
        loader.iterate(collector)
        return collector.snps, collector.snps_by_rsid
    
class ProcessWeightDB(object):
    def __init__(self, args):
        # input and output 
        self.utmost_dir = args.utmost_dir
        self.weight_db = pathLeaf(args.weight_db)
        self.db_path = args.weight_db
        self.data_folder = args.input_folder
        self.gene_list = args.gene_list.split(",")
        test_str = "_".join(self.gene_list) + "_" + args.gwas_str        
        self.covariance_output = self.utmost_dir + 'intermediate/test/' + test_str + "/"
        self.input_format = Formats.PrediXcan
        self.db_logic_gene = {}
        self.db_logic_tissue = {}
        self.found_genes_for_covariance = {}
        self.found_genes_for_correlation = {}
        self.db_logic_dict = {}
        self.gene_rsid_dict = {}
        self.db_file_list = []
        self.gene_count = 0
        self.chr_idx = args.chr_idx
        self.gwas_df = None
        
        # preprocessing filter
        self.min_maf_filter = float(args.min_maf_filter) if args.min_maf_filter else None
        self.max_maf_filter = float(args.max_maf_filter) if args.max_maf_filter else None
        self.max_snps_in_gene = int(args.max_snps_in_gene) if args.max_snps_in_gene else None

    def run(self):
        # run the main function
        if not self.covariance_output:
            logging.info("Provide --covariance_output or both")
            return

        
        # list all the databases in the path
        for file in sorted(os.listdir(self.db_path)):
            if file.endswith(".db"):
                self.db_file_list.append(file)
                
        
        # load the database and build the separate db entry logic  
        logging.info("Loading Weights")
        count = 0    
        # gene level snplist
        for gene in self.gene_list:
            self.gene_rsid_dict[gene] = []
        # reading weight database
        tissue_vec = []
        #self.db_file_list = self.db_file_list[0:1]
        for file in self.db_file_list:
            count += 1
            filename = self.db_path + file
            tissue = file.split(".db")[0]
            tissue_vec.append(tissue)
            self.db_logic_dict[tissue] = WeightDBUtilities.WeightDBEntryLogicGene(filename, self.gene_list)
            self.db_logic_tissue[tissue] = copy.deepcopy(self.db_logic_dict[tissue])
            logging.info("Building file" + str(count))
            print "INFO: Building file" + str(count)
            #print self.db_logic_dict[tissue].weights_by_gene[self.gene_list[0]]
        
        # same tissue, different genes
        for tissue in tissue_vec:
            tmp_db_logic = self.db_logic_tissue[tissue]
            tmp_db_logic.weights_by_gene['merged'] = {}
            for gene in self.gene_list:
                if gene in tmp_db_logic.weights_by_gene.keys():
                    tmp_db_logic.weights_by_gene['merged'].update(tmp_db_logic.weights_by_gene[gene])
                    del tmp_db_logic.weights_by_gene[gene]
            if len(tmp_db_logic.weights_by_gene['merged'].keys()) > 0:
                rsid_merged = tmp_db_logic.weights_by_gene['merged'].keys()
                for rsid in rsid_merged:
                    tmp_db_logic.genes_for_an_rsid[rsid] = 'merged'
                # build related data
                self.buildFilesTissue(tmp_db_logic, tissue, self.chr_idx)
                
        # same gene, different tissues
        for gene in self.gene_list:
            self.db_logic_gene[gene] = copy.deepcopy(self.db_logic_dict[tissue_vec[0]])
            tmp_db_logic = self.db_logic_gene[gene]
            if not gene in tmp_db_logic.weights_by_gene.keys():
                tmp_db_logic.weights_by_gene[gene] = {}
            for tissue in tissue_vec:
                if gene in self.db_logic_dict[tissue].weights_by_gene.keys():
                    tmp_db_logic.weights_by_gene[gene].update(self.db_logic_dict[tissue].weights_by_gene[gene])
            for tmp_gene in self.gene_list:
                if tmp_gene != gene and tmp_gene in tmp_db_logic.weights_by_gene.keys():
                    del tmp_db_logic.weights_by_gene[tmp_gene]
            if len(tmp_db_logic.weights_by_gene[gene].keys()) > 0:
                rsid_merged = tmp_db_logic.weights_by_gene[gene].keys()
                for rsid in rsid_merged:
                    tmp_db_logic.genes_for_an_rsid[rsid] = gene  
                # build related data
                self.buildFilesGene(tmp_db_logic, gene, self.chr_idx) 
            else:
                print "Not coresponding SNPs for gene " + gene
                return      
       
        # summary of gene count and snp count
        logging.info("Total Genes:" + str(len(self.db_logic_dict[tissue].weights_by_gene.keys())))

        logging.info("Preprocess successfully")
        print "INFO: Preprocess complete"

    def passGWAS(self, gwas_df):
        self.gwas_df = gwas_df
        
    def checkGWASFlipping(self, ):
        self.gwas_df = gwas_df

    def getSNPS(self, name, weight_db_logic):
        dosageLoader = None
        if self.input_format == Formats.IMPUTE:
            dosageLoader = ThousandGenomesUtilities.IMPUTEDosageLoader(self.data_folder, name) #outdated code
        elif self.input_format == Formats.PrediXcan:
            dosageName = Utilities.dosageName(name)
            path = os.path.join(self.data_folder, dosageName)
            dosageLoader = PrediXcanFormatDosageLoader(path, weight_db_logic)
        else:
            logging.info("Invalid input format: %s", self.input_format)
            return
        snps, snps_by_rsid = dosageLoader.load()
        #print snps
        return snps, snps_by_rsid
    
    # build the covariance file for each tissue
    def buildFilesGene(self, weight_db_logic, gene, chr_idx):
        weight_db_logic = self.db_logic_gene[gene]
        do_covariances = self.covariance_output is not None
        if do_covariances:
            covariance_dir = os.path.dirname(self.covariance_output)
            if not os.path.exists(covariance_dir):
                os.makedirs(covariance_dir)

        if not do_covariances:
            return
        # get snp info
        name = 'EUR.chr' + str(chr_idx)
        snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)
        #print "calculating cov..."
        self.addToCovarianceFile(weight_db_logic, name, snps, snps_by_rsid)

    def addToCovarianceFile(self, weight_db_logic, name, snps, snps_by_rsid):
        logging.info("Adding to covariance for %s-%s", name, self.weight_db)
        # get the total genes in that weight_db_logic
        genes = weight_db_logic.weights_by_gene.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            print(gene)
            if percent == last_reported_percent+1:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent
            self.buildCovarianceEntries(name, gene, weight_db_logic, snps_by_rsid)
    
    def buildCovarianceEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        # get a dict of specific gene
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = weights_in_gene.keys()
        #gather as much data as we can work on
        related_rsids, related_data, flip_tag_vec = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)
        if len(related_rsids) == 0:
            return []
        #covariance matrix of related SNP's data
        array = np.array(related_data)
        cov = np.cov(array)
        self.gene_count += 1;
        logging.info("GENE:" + str(self.gene_count))
        print "INFO: build covariance"
        # save the covariance matrix
        # write the cov matrix
        cov_dir = self.covariance_output 
        cov_filename = cov_dir + gene + ".cov"
        if not os.path.exists(cov_dir):
            os.mkdir(cov_dir)
        with open(cov_filename,"w") as fo:
            if cov.shape == ():
                fo.write(str(cov) + "\n")
            else:  
                np.savetxt(fo, cov)

        snp_filename = self.covariance_output + gene + ".snplist"
        with open(snp_filename,"w") as fo:
            for snp in related_rsids:
                fo.write(snp + "\n")
                
        # build gwas
        gwas_output = self.covariance_output + '/' + gene + '.gene.gwas'
        buildGWAS(gwas_df, related_rsids, gwas_output, flip_tag_vec)
        
    def buildRelatedData(self, rsids_from_genes, snps_by_rsid, weights_in_gene):
        related_rsids = []
        related_data = []
        flip_tag_vec = []

        l = len(rsids_from_genes)
        if self.max_snps_in_gene and l > self.max_snps_in_gene:
            logging.info("Skipping covariance too large: %d", l)
            return related_data, related_rsids, flip_tag_vec

        for rsid in rsids_from_genes:
            if not rsid in snps_by_rsid:
                print "related rsid not present in genotype data"
                logging.log(5, "related rsid %s not present in genotype data", rsid)
                continue

            related_snp = snps_by_rsid[rsid]
            freq = sum(related_snp.data)*1.0/(2*len(related_snp.data))
            if self.min_maf_filter and self.min_maf_filter > freq:
                logging.log(6, "related rsid %s below min maf: %s", rsid, freq)
                continue

            if self.max_maf_filter and self.max_maf_filter < freq:
                logging.log(6, "related rsid %s  above max maf: %s", rsid, freq)
                continue
            # dosage data
            data = related_snp.data
            weight = weights_in_gene[rsid]
            if weight.ref_allele == related_snp.eff_allele and\
            weight.eff_allele == related_snp.ref_allele:
                logging.log(7, "related rsid %s has alleles flipped compared to model, transforming dosage", rsid)
                data = map(lambda x: 2-x, data)
            related_data.append(data)
            related_rsids.append(rsid)
            
        # GWAS flipping
        tmp_df = copy.deepcopy(self.gwas_df[self.gwas_df['snp'].isin(related_rsids)])
        for index, row in tmp_df.iterrows():
            flip_tag = 1
            rsid = row['snp']
            eff_allele = row['effect_allele']
            ref_allele = row['non_effect_allele']
            weight = weights_in_gene[rsid]
            if weight.ref_allele == eff_allele and weight.eff_allele == ref_allele:
                flip_tag = -1
            flip_tag_vec.append(flip_tag)
        return related_rsids, related_data, flip_tag_vec
        
     
    # build the covariance file for each tissue
    def buildFilesTissue(self, weight_db_logic, tissue, chr_idx):
        weight_db_logic = self.db_logic_tissue[tissue]
        do_covariances = self.covariance_output is not None
        if do_covariances:
            covariance_dir = os.path.dirname(self.covariance_output)
            if not os.path.exists(covariance_dir):
                os.makedirs(covariance_dir)

        if not do_covariances:
            return
        # get snp info
        name = 'EUR.chr' + str(chr_idx)
        snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)
        #print "calculating cov..."
        self.addToCovarianceFileTissue(weight_db_logic, name, snps, snps_by_rsid, tissue)
        
    def addToCovarianceFileTissue(self, weight_db_logic, name, snps, snps_by_rsid, tissue):
        logging.info("Adding to covariance for %s-%s", name, self.weight_db)
        # get the total genes in that weight_db_logic
        genes = weight_db_logic.weights_by_gene.keys()
        total_genes = len(genes)
        last_reported_percent = 0
        processed = 0
        for gene in genes:
            processed += 1
            percent = int(processed*100.0 / total_genes)
            print "INFO: " + gene
            if percent == last_reported_percent+1:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent
            self.buildCovarianceEntriesTissue(name, gene, weight_db_logic, snps_by_rsid, tissue)
            
    def buildCovarianceEntriesTissue(self, name, gene, weight_db_logic, snps_by_rsid, tissue):
        # get a dict of specific gene
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = weights_in_gene.keys()
        # gather as much data as we can work on
        related_rsids, related_data, flip_tag_vec = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)

        if len(related_rsids) == 0:
            return []

        #covariance matrix of related SNP's data
        array = np.array(related_data)
        cov = np.cov(array)
        self.gene_count += 1;
        logging.info("GENE:" + str(self.gene_count))
        print "INFO: build covariance"
        # save the covariance matrix
        # write the cov matrix
        cov_dir = self.covariance_output + "/" + tissue + "/" 
        cov_filename = cov_dir + gene + ".cov"
        if not os.path.exists(cov_dir):
            os.mkdir(cov_dir)
        with open(cov_filename,"w") as fo:
            if cov.shape == ():
                fo.write(str(cov) + "\n")
            else:  
                np.savetxt(fo, cov)

        snp_filename = self.covariance_output + "/" + tissue + "/" + gene + ".snplist"
        with open(snp_filename,"w") as fo:
            for snp in related_rsids:
                fo.write(snp + "\n")
        
        # build gwas
        gwas_output = self.covariance_output + '/' + tissue + '/' + 'merged.gwas'
        #print flip_tag_vec
        buildGWAS(gwas_df, related_rsids, gwas_output, flip_tag_vec)

class conditionalTest(object):
    def __init__(self, args):
        # input and output 
        self.utmost_dir = args.utmost_dir
        self.weight_db = args.weight_db
        self.db_path = args.weight_db
        self.data_folder = args.input_folder
        self.gene_list = args.gene_list.split(",")
        test_str = "_".join(self.gene_list) + "_" + args.gwas_str        
        self.covariance_output = self.utmost_dir + 'intermediate/test/' + test_str + "/"
        self.input_format = Formats.PrediXcan
        self.outcome = None
        self.zscore_gene = None
        self.causal_dir = self.covariance_output
        self.tmp_dir = self.covariance_output + "/temp/"
        self.output_dir = args.output_dir
        self.tissue_list = []

    def run(self):
        fi = []
        for file in sorted(os.listdir(self.weight_db)):
            if file.endswith(".db"):
                fi.append(file)
        
        gene_list = self.gene_list
        n = len(fi)
        p = len(self.gene_list)
        
        # initialize outcome matrix
        #outcome = np.zeros(shape=(p, n+1))
        #zscore_gene = np.zeros(shape=(p, n))
        self.outcome = np.full((p, n+1), np.nan)
        self.zscore_gene = np.full((p, n), np.nan)
        outcome = self.outcome
        zscore_gene = self.zscore_gene
        
        # build tissue weights dict
        gene_tissue_weights = {}
        for cur_gene in gene_list:
            gene_dir = self.covariance_output + "/" + cur_gene + ".snplist"
            gene_rsid = pd.read_csv(gene_dir, sep = " ", header = None)[0].tolist()
            gene_tissue_weights[cur_gene] = np.empty((0,len(gene_rsid)))
        
        # read in preprocessed data
        print "INFO: Reading preprocessed data"
        for tissue_idx, dbname in enumerate(fi):
            tissue_name = dbname.split(".")[0]
            self.tissue_list.append(tissue_name)
            # connect to the database
            conn = create_connection(self.weight_db + dbname)
            cur = conn.cursor()  
            print "INFO: "  + str(round(float(tissue_idx + 1)/float(n) * 100, 2)) + "%"
            # corresponding paths
            gene_tissue_str = self.causal_dir + "/" + tissue_name
            #gene_level_str = self.causal_dir
            snp_dir = gene_tissue_str + "/merged.snplist"
            cov_dir = gene_tissue_str + "/merged.cov"
            gwas_dir = gene_tissue_str + "/merged.gwas"    
            if not os.path.exists(gwas_dir):
                continue
            gwas_df = pd.read_csv(gwas_dir, sep = "\t")
            # read snp list data
            rsid_df = pd.read_csv(snp_dir, sep = " ", header = None)
            snp_rsid = list(rsid_df.loc[:,0])  
            # read cov data
            cov_snp = pd.read_csv(cov_dir, sep = " ", header = None)    
            # matrix of weights
            m = len(snp_rsid)
            weights = np.zeros(shape = (m, p)) 
            for i in range(p):
                sql_q = 'select * from weights where gene = "' + gene_list[i] + '"'
                # return all the items
                tmp_query = cur.execute(sql_q).fetchall()
                # extract the snp
                rsid_in_db = list(map(lambda x: str(x[0]), tmp_query))
                match_idx = match_list(rsid_in_db, snp_rsid)
                indi = match_idx[match_idx > -1]
                # extract the wegihts
                #tmp_weights = np.array(list(map(lambda x: str(x[2]), tmp_query)))
                tmp_weights = np.array(list(map(lambda x: (x[2]), tmp_query)))
                if sum(match_idx > -1) > 0:
                    weights[indi, i] = tmp_weights[match_idx > -1]  
            # covariance matrix of genes
            cov_genes = np.zeros(shape = (p, p))
            for a in range(p):
                for b in range(p):
                    cov_genes[a, b] = (np.transpose(weights[:, a])).dot(cov_snp).dot(weights[:, b])
            #print(cov_genes)
            gene_appear = np.where(np.diag(cov_genes) != 0)[0]
            #print(gene_appear[0])
            if len(gene_appear) == 0:
                continue
            #gene_appear = np.array(gene_appear)
            cov_gene_inv = np.linalg.inv(cov_genes[gene_appear, :][:, gene_appear])
            sd = np.sqrt(np.diag(cov_snp))   
            # exatract zscore from gwas file
            zscore = np.zeros(m)
            match_idx = match_list(snp_rsid, list(gwas_df.iloc[:, 0]))
            indi = match_idx[match_idx > -1]
            zscore[match_idx > -1] = gwas_df["zscore"].iloc[indi].values 
            #print(zscore)
            # zscore gene
            if len(sd) == 1:
                gene_weights =cov_gene_inv.dot(np.transpose(weights[:, gene_appear])).dot(sd)/np.sqrt(np.diag(cov_gene_inv))
                zscore_gene[gene_appear, tissue_idx] = cov_gene_inv.dot(np.transpose(weights[:, gene_appear])).dot(sd) \
                .dot(zscore)/np.sqrt(np.diag(cov_gene_inv))
            else:
                gene_weights =cov_gene_inv.dot(np.transpose(weights[:, gene_appear])).dot(np.diag(sd))/(np.sqrt(np.diag(cov_gene_inv))[:,None])
                zscore_gene[gene_appear, tissue_idx] = cov_gene_inv.dot(np.transpose(weights[:, gene_appear])).dot(np.diag(sd)) \
                .dot(zscore)/np.sqrt(np.diag(cov_gene_inv))
            for gene_idx in range(len(gene_appear)):
                cur_gene = gene_list[int(gene_appear[gene_idx])]
                gene_dir = self.covariance_output + "/" + cur_gene + ".snplist"
                gene_rsid = pd.read_csv(gene_dir, sep = " ", header = None)[0].tolist()
                sql_q = 'select * from weights where gene = "' + cur_gene + '"'
                # return all the items
                tmp_query = cur.execute(sql_q).fetchall()
                # extract the snp
                gene_tissue_rsid = list(map(lambda x: str(x[0]), tmp_query))
                index_1 = match_list(gene_tissue_rsid, snp_rsid)
                index_2 = match_list(gene_tissue_rsid, gene_rsid)
                indi_1 = index_1 > -1
                indi_2 = index_2 > -1
                tmp_gene_tissue_weights = np.zeros(len(gene_rsid))
                #output_file = self.tmp_dir + cur_gene + ".test.gene.weights"
                if len(index_2) != 0:
                    if gene_weights.shape == (1,):
                        gene_weights = gene_weights.reshape(1, -1)
                    tmp_gene_tissue_weights[index_2[indi_2]] = gene_weights[gene_idx, index_1[indi_1]]
                    gene_tissue_weights[cur_gene] = np.vstack((gene_tissue_weights[cur_gene],tmp_gene_tissue_weights))
                outcome[gene_appear, tissue_idx] = 2 * st.norm.cdf(-np.absolute(zscore_gene[gene_appear, tissue_idx]))
        
        # output the gene tissue weights
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        for cur_gene in gene_list:
            np.savetxt(self.tmp_dir + cur_gene + ".test.gene.weights", gene_tissue_weights[cur_gene])
        
        # GBJ test
        print "INFO: GBJ test"
        for gene_idx, cur_gene in enumerate(gene_list):
            index = ~np.isnan(outcome[gene_idx, 0:n])
            if (sum(index) > 1):
                output_file = self.tmp_dir + cur_gene + ".test.gene.weights"
                #output_file_w = causal_dir + tar_gene + "/" + cur_gene + ".test.gene.weights"
                weights_gene = pd.read_csv(output_file, sep = " ", header = None)
                output_file_cov = self.covariance_output + "/" + cur_gene + ".cov"
                #print(weights_gene)
                #print(cov_snp)

                cov_snp = pd.read_csv(output_file_cov, sep = " ", header = None)
                cov_gene = weights_gene.dot(cov_snp).dot(np.transpose(weights_gene))

                # normalization
                for i in range(cov_gene.shape[0]):
                    if cov_gene.iloc[i, i] != 0:
                        cov_gene.iloc[i, :] = cov_gene.iloc[i, :]/np.sqrt(cov_gene.iloc[i, i])
                        cov_gene.iloc[:, i] = cov_gene.iloc[:, i]/cov_gene.iloc[i, i]
                #print(cov_gene)
                if np.allclose(cov_gene,cov_gene.T):
                    # GBJ
                    cov_gene = cov_gene.round(8)
                    np.savetxt(self.tmp_dir + cur_gene + ".cov", cov_gene)
                    np.savetxt(self.tmp_dir + cur_gene + ".zscore", zscore_gene[gene_idx, index])
                    return_val = subprocess.call(["Rscript",  self.utmost_dir + "gbj.R", cur_gene, self.tmp_dir])

                    # run the test            
                    #p_value = np.read_csv(tmp_dir + cur_gene + ".pvalue", header = None)
                    with open(self.tmp_dir + cur_gene + ".pvalue", 'r') as fi:
                        p_value = fi.readline().strip()
                    # output the test result to the result matrix
                    outcome[gene_idx, n] = p_value
                    #print(GBJ_res.rx2("GBJ_pvalue")[0])
            elif sum(index) == 1:
                outcome[gene_idx, n] = outcome[gene_idx, 0:n][index]
            else:
                outcome[gene_idx, n] = np.nan
                    
        # output the results
        self.tissue_list.append("Overall")
        output_df = pd.DataFrame(outcome, index = self.gene_list, columns= self.tissue_list)
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        output_df.to_csv(self.output_dir + "results.csv", na_rep = "NA")
        #np.savetxt(self.output_dir + "results.txt", outcome)
        print "INFO: Conditional test complete!" 
        print "INFO: " + self.output_dir + "results.csv"
   

if __name__ == "__main__":
    # argument processing
    import argparse
    parser = argparse.ArgumentParser(description='Perform cross tissue conditional analysis.')

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--utmost_dir",
                        help="path of utmost",
                        default="./")
    
    parser.add_argument("--weight_db",
                        help="name of weight db in data folder",
                        default="sample_data/weigth_db_GTEx/")
    
    parser.add_argument("--covariance_output",
                        help="the covariance directory",
                        default="intermediate/")
    
    parser.add_argument("--gene_list",
                        help="input gene list",
                        default=None)

    parser.add_argument("--input_folder",
                        help="name of folder containing dosage data",
                        default="sample_data/dosage/")
    
    parser.add_argument("--output_dir",
                        help="the name of output file",
                        default="./outcome/")
    
    parser.add_argument("--gwas_str",
                        help="name of the gwas study",
                        default="gwas")
    
    parser.add_argument("--gwas_folder",
                        help="name of folder containing GWAS data. All files in the folder are assumed to belong to a single study.",
                        default="sample_data/GWAS")

    parser.add_argument("--gwas_file_pattern",
                        help="Pattern to recognice GWAS files in folders (in case there are extra files and you don't want them selected).",
                        default=None)

    GWASUtilities.add_gwas_arguments_to_parser(parser)

    parser.add_argument("--separator",
                        help="Character or string separating fields in input file. Defaults to any whitespace.",
                        default=None)
    
    parser.add_argument("--skip_until_header",
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.",
                        default=None)
    
    parser.add_argument('--min_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument('--max_maf_filter',
                   help="Filter snps according to this maf",
                   default=None)

    parser.add_argument("--max_snps_in_gene",
                        help="Ignore any gene that has snps above this value",
                        type=int,
                        default=None)
    
    parser.add_argument("--chr_idx",
                        help="Chromosome index of input genes",
                        default=None)
                        
    
    args = parser.parse_args()
    
    # reading GWAS
    gwas_df = readGWAS(args)

    # data preprocessing
    pwdb = ProcessWeightDB(args)
    pwdb.passGWAS(gwas_df)
    pwdb.run()

    # conditional test
    ct = conditionalTest(args)
    ct.run()
