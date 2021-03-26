library(MASS)
library(spatstat)
library(cfdr)
library(ggplot2)
library(cowplot)
library(fields)
library(dplyr)
library(polyCub)
library(hexbin)
library(bigsplines)
library(locfdr)
library(data.table)
library(fcfdr)

df <- readRDS("asthma_p.RDS") # data.frame of rsIDs (SNPID) and asthma GWAS p-values (european_ancestry_pval_rand)

# independent subset of SNPs
ind_snps <- readRDS("ind_snps.RDS") # names of SNPs in the independent subset of SNPs
indep_index <- na.omit(match(ind_snps, df$SNPID))
orig_p <- df$european_ancestry_pval_rand

q_vals <- readRDS("../H3K27ac_data/H3K27ac_qs.RDS")

# make containers for V values (with first column original P values for first iteration)
v_vals <- data.frame("p" = orig_p, "v1" = orig_p, "v2" = orig_p)

for(i in 1:2){
  
  p <- v_vals[,i]
  q <- q_vals[,i]
  
  print(paste0("Starting iteration ",i," p is ", p[1]," q is ", q[1]))
  
  # run flexible cFDR
  res_full <- flexible_cfdr(p = p, q = q, indep_index, splinecorr = FALSE)
  
  saveRDS(res_full, paste0("res_iter",i,".RDS"))
  
  res <- res_full[[1]]
  v_vals[,i+1] <- res$v
  q_vals[,i] <- res$q # replace to match sign used in func cFDR
  
  print(paste0("Finished iteration ", i, " v is ", res$v[1]))
  
  rm(res)
  
}

saveRDS(data.frame(v_vals, q_vals), "H3K27ac_cFDR_res.RDS")

### Boca and Leek's FDR regression

library("swfdr")

GWAS_lm_qvalue <- lm_qvalue(orig_p, X = q_vals[,c(1,2)])

saveRDS(GWAS_lm_qvalue, "H3K27ac_BL_res.RDS")
