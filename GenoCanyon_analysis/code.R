library(MASS)
library(spatstat)
library(cfdr)
library(fields)
library(ggplot2)
library(cowplot)
library(dplyr)
library(polyCub)
library(hexbin)
library(bigsplines)
library(locfdr)
library(fcfdr)

df <- readRDS("asthma_p.RDS") # data.frame of rsID (SNPID) and asthma GWAS p-values (european_ancestry_pval_rand) for 1968651 SNPs in asthma analysis
ind_snps <- readRDS("ind_snps.RDS") # vector of rsIDs for the independent subset of SNPs
indep_index <- na.omit(match(ind_snps, df$SNPID))
orig_p <- df$european_ancestry_pval_rand

q_tmp <- readRDS("GenoCanyon_Prediction.RDS") # data.frame of chr (V1), pos (V2) and GenoCanyon score (V3) for 1968651 SNPs in asthma analysis
q <- q_tmp$V3

# run flexible cFDR
res <- flexible_cfdr(p = orig_p, q, indep_index)

saveRDS(res, "GC_func_cFDR_res.RDS")