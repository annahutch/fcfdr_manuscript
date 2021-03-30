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
library(IHW)
library(swfdr)

df <- readRDS("../asthma_p.RDS") # data.frame of rsID (SNPID) and asthma GWAS p-values (european_ancestry_pval_rand) for 1968651 SNPs in asthma analysis
ind_snps <- readRDS("../ind_snps.RDS") # vector of rsIDs for the independent subset of SNPs
indep_index <- na.omit(match(ind_snps, df$SNPID))
orig_p <- df$european_ancestry_pval_rand

q_tmp <- readRDS("GenoCanyon_Prediction.RDS") # data.frame of chr (V1), pos (V2) and GenoCanyon score (V3) for 1968651 SNPs in asthma analysis
q <- q_tmp$V3

fdr_thr <- 0.000148249 # FDR threshold corresponding to conventional genome-wide significance p-value threshold

# run IHW
ihw_res <- ihw(orig_p, q, alpha = fdr_thr)

# run Boca and Leek's FDR regression
BL_res <- lm_qvalue(orig_p, X = q)

# run flexible cFDR
cfdr_res <- flexible_cfdr(p = orig_p, q, indep_index, maf = df$MAF)

saveRDS(ihw_res, "GC_IHW_res.RDS")
saveRDS(BL_res, "GC_BL_res.RDS")
saveRDS(cfdr_res, "GC_func_cFDR_res.RDS")

