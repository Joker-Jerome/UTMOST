library(GBJ)

# reading 
commandArgs(trailingOnly = TRUE)
args <- commandArgs(TRUE)
gene <- args[1] 
dir <- args[2] 

# cov zscore
cov <- read.table(paste0(dir, gene, ".cov"), header = F)
zscore <- read.table(paste0(dir, gene, ".zscore"), header = F)
cov <- data.matrix(cov)
zscore <- zscore$V1
dimnames(cov) <- NULL

#GBJ test and results
res <- GBJ(test_stats=zscore, cor_mat=cov)
write.table(res$GBJ, paste0(dir, gene, ".teststat"), row.names = F, col.names = F, quote = F)
write.table(res$GBJ_pvalue, paste0(dir, gene, ".pvalue"), row.names = F, col.names = F, quote = F)


