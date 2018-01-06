#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Covariance Calculation-- 


It defines classes_and_methods

@author:     heroico, zhao_lab 

@copyright:  2017 organization_name. All rights reserved.

@license:    license

@contact:    zhao_lab at Yale
@deffield    updated: Updated
'''

import logging
import numpy as np
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


def pathLeaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def mergeTwoDicts(x, y):
    z = x.copy()  
    z.update(y)   
    return z

class ProcessWeightDB(object):
    def __init__(self, args):
        # input and output 
        self.weight_db = pathLeaf(args.weight_db)
        self.db_path = args.weight_db
        self.data_folder = args.input_folder
        self.covariance_output = args.covariance_output
        if args.covariance_output is None:
            comp = os.path.splitext(self.weight_db)[0]
            name = comp + ".cov.txt.gz"
            path = os.path.join("intermediate", "cov")
            path = os.path.join(path, name)
            self.covariance_output = path    
        self.input_format = args.input_format
        self.store_pickle_only = args.store_pickle_only
        
        # data entry logic
        self.found_genes_for_covariance = {}
        self.found_genes_for_correlation = {}
        self.db_logic_dict = {}
        self.gene_rsid_dict = {}
        self.db_file_list = []
        
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
            if file.endswith(".db") and not file.endswith("sqtl.db"):
                self.db_file_list.append(file)
        
        # load the database and build the separate db entry logic
        logging.info("Loading Weights")
        count = 0      
        for file in self.db_file_list:
            count += 1
            filename = self.db_path + file
            self.db_logic_dict[file] = WeightDBUtilities.WeightDBEntryLogic(filename)
            logging.info("Building file" + str(count))
        
        # merge the info from different databases
        tmp_logic_object = self.db_logic_dict[list(self.db_logic_dict.keys())[0]]
        count = 0
        for db_logic in self.db_logic_dict.values():
            count += 1
            logging.info("Scanning file" + str(count))
            # update the weights_by_gene  
            count_gene = 0
            num_gene = len(db_logic.weights_by_gene.keys())
            for gene in db_logic.weights_by_gene.keys():
                count_gene += 1 
                if count_gene % 150 == 0:
                    logging.info("Percentage of genes processed  " + str(round(float(count_gene)/float(num_gene),2)*100))
                if gene in tmp_logic_object.weights_by_gene.keys():
                    for rsid in db_logic.weights_by_gene[gene].keys():
                        if rsid not in tmp_logic_object.weights_by_gene[gene].keys():
                            tmp_logic_object.weights_by_gene[gene][rsid] = db_logic.weights_by_gene[gene][rsid]              
                else:
                    tmp_logic_object.weights_by_gene[gene] = db_logic.weights_by_gene[gene]             
                
        # summary of gene count and snp count
        logging.info("Total Genes:" + str(len(tmp_logic_object.weights_by_gene.keys())))
        rsid_count = 0
        for gene in tmp_logic_object.weights_by_gene.keys():
            rsid_count += len(tmp_logic_object.weights_by_gene[gene].keys())
        logging.info("Total SNPs:" + str(rsid_count))
    
        # store the pickle file
        pickle_out = open("db_weight_logic.pickle","wb")
        pickle.dump(tmp_logic_object,pickle_out)
        pickle_out.close()   
        
        # store the gene info
        self.saveGeneInfo(tmp_logic_object)
            
        # whether calculate the covariance directly
        # store the database entry logic as pickle file
        if not self.store_pickle_only:
            self.buildFiles(tmp_logic_object)

        logging.info("Ran successfully")
    
    # build the covariance file
    def buildFiles(self, weight_db_logic):

        do_covariances = self.covariance_output is not None
        if do_covariances:
            if os.path.exists(self.covariance_output):
                logging.info("%s already exists, delete it if you want it figured out again", self.covariance_output)
                do_covariances = False
            else:
                covariance_dir = os.path.dirname(self.covariance_output)
                if not os.path.exists(covariance_dir):
                    os.makedirs(covariance_dir)
                self.writeFileHeader(self.covariance_output)

        if not do_covariances:
            return
        # load the dosage data from dosage folder
        names = Utilities.dosageNamesFromFolder(self.data_folder)
        for name in names:
            snps, snps_by_rsid = self.getSNPS(name, weight_db_logic)
            if do_covariances:
                self.addToCovarianceFile(weight_db_logic, name, snps, snps_by_rsid)

    def getSNPS(self, name, weight_db_logic):
        dosageLoader = None
        if self.input_format == Formats.IMPUTE:
            dosageLoader = ThousandGenomesUtilities.IMPUTEDosageLoader(self.data_folder, name) #outdated code
        elif self.input_format == Formats.PrediXcan:
            dosageName = Utilities.dosageName(name)
            path = os.path.join(self.data_folder, dosageName)
            dosageLoader = PrediXcanFormatUtilities.PrediXcanFormatDosageLoader(path, weight_db_logic)
        else:
            logging.info("Invalid input format: %s", self.input_format)
            return
        snps, snps_by_rsid = dosageLoader.load()
        return snps, snps_by_rsid

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
            print(len(weight_db_logic.weights_by_gene[gene].keys()))
            if percent == last_reported_percent+1:
                logging.info("%d percent genes processed", percent)
                last_reported_percent = percent

            self.buildCovarianceEntries(name, gene, weight_db_logic, snps_by_rsid)

    def buildCovarianceEntries(self, name, gene, weight_db_logic, snps_by_rsid):
        # get a dict of specific gene
        weights_in_gene = weight_db_logic.weights_by_gene[gene]
        rsids_from_genes = weights_in_gene.keys()
        
        #gather as much data as we can work on
        related_rsids, related_data = self.buildRelatedData(rsids_from_genes, snps_by_rsid, weights_in_gene)

        if len(related_rsids) == 0:
            return []

        self.updateFoundCovariance(gene, name)

        #covariance matrix of related SNP's data
        array = np.array(related_data)
        cov = np.cov(array)
        self.total_gene += 1;
        logging.info("GENE:" + str(self.total_gene))
        
        #save the covariance matrix
        # write the cov matrix
        cov_filename = self.covariance_output + "/" + gene + ".cov"
        with open(cov_filename,"w") as fo:
            if cov.shape == ():
                fo.write(str(cov) + "\n")
            else:  
                np.savetxt(fo, cov)
            
        snp_filename = self.covariance_output + "/" + gene + ".snplist"
        with open(snp_filename,"w") as fo:
            for snp in related_rsids:
                fo.write(snp + "\n")

    def buildRelatedData(self, rsids_from_genes, snps_by_rsid, weights_in_gene):
        related_rsids = []
        related_data = []

        l = len(rsids_from_genes)
        if self.max_snps_in_gene and l > self.max_snps_in_gene:
            logging.info("Skipping covariance too large: %d", l)
            return related_data, related_rsids

        for rsid in rsids_from_genes:
            if not rsid in snps_by_rsid:
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
        return related_rsids, related_data
    
    # build the covariance file
    def saveGeneInfo(self, weight_db_logic):
        gene_list = list(weight_db_logic.weights_by_gene.keys())
        with open("gene_info.txt","w") as fo:
            fo.write("gene_ensg" + "\n")
            for item in gene_list:
                fo.write(item + "\n")       


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Build covariances from dosage data and weights database.')

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    parser.add_argument("--weight_db",
                        help="name of weight db in data folder",
                        default="data/DGN-WB_0.5.db")

    parser.add_argument("--input_folder",
                        help="name of folder containing dosage data",
                        default="intermediate/TGF_EUR")

    parser.add_argument("--covariance_output",
                        help="Name of file to dump covariance results in. Defaults to './intermediate/cov/' + file name prefix from '--weight_db' argument",
                        default="./intermediate/cov/")

    parser.add_argument('--input_format',
                   help='Input dosage files format. Valid options are: IMPUTE, PrediXcan',
                   default=Formats.PrediXcan)
    

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
    
    parser.add_argument("--store_pickle_only",
                        help="Only store the database entry logic",
                        type=int,
                        default=0)

    args = parser.parse_args()

    Logging.configureLogging(int(args.verbosity))
    
    
    work = ProcessWeightDB(args)
    work.run()  
    
    

