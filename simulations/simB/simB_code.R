##### code for related trait simulations (simulation B)

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

# simulated GWAS data
files <- list.files(path = "../simulatedGWAS", pattern="*.RDS", full=TRUE)

# LDAK weights for all SNPs  
# to find independent subset of SNPs
sim_weights <- readRDS("../sim_weights.RDS")
indep_index <- which(sim_weights != 0)

func <- function(){
  
  data <- lapply(files, function(x){ # randomly select a simulations from each LD block
    data_full <- readRDS(x)
    tmp <- data_full[[sample(1:length(data_full), 1)]]
    subset(tmp, select = c("z","p","CV","pos","block","max_r2"))
  }) %>% rbindlist()
  
  p <- data$p
  block <- data$block
  
  # make containers for V values (with first column original P values for first iteration)
  v_vals <- data.frame("p" = p, "v1" = p, "v2" = p, "v3" = p, "v4" = p, "v5" = p)
  v_vals_emp <- data.frame("p_emp" = p, "v1_emp" = p, "v2_emp" = p, "v3_emp" = p, "v4_emp" = p, "v5_emp" = p)
  
  # simulate q as p-vals from related traits
  # holder
  q_vals <- data.frame("q1" = runif(length(p), 0, 1), "q2" = runif(length(p), 0, 1), "q3" = runif(length(p), 0, 1), "q4" = runif(length(p), 0, 1), "q5" = runif(length(p), 0, 1))
  
  niter <- 5
  
  # make empty list for shared CVs in each iter
  shared_CVs <- vector(mode = "list", length = 5)
  
  for(j in 1:niter){
    
    set.seed(j*1000)
    
    p <- v_vals[,j]
    
    q <- q_vals[,j]
    
    shared_CV <- 0
    
    while(length(shared_CV)<20 & cor(p,q)<0.15){ # ensure correlation is high enough between "related" traits
      
      q_out_tmp <- lapply(files, function(x){ # randomly select an element from each block
        data_full <- readRDS(x)
        subset <- data_full[[sample(1:length(data_full), 1)]]
        data.frame(q = subset$p, CV = subset$CV)
      }) %>% rbindlist()
      
      q <- q_out_tmp[[1]]
      q_CV <- q_out_tmp[[2]]
      
      # shared CVs
      shared_CV <- intersect(which(q_CV==1), which(data$CV==1))
      
    }
    
    q_vals[,j] <- q
    
    print(paste0("Iteration: ", j))
    print(paste0("Correlation is: ",cor(p, q)))
    print(paste0("Ind correlation is: ",cor(p[indep_index], q[indep_index])))
    print(paste0("q is :", c(q[1])))
    print(paste0("no shared CVs is :", length(shared_CV)))
    
    # run flexible cFDR
    res_full <- flexible_cfdr(p = p, q = q, indep_index, gridp = 0) # no left censoring required due to distribution of aux data
    
    res <- res_full[[1]]
    v_vals[,j+1] <- res$v
    
    # run empirical cFDR with fold removal
    fold <- vector(mode = "list", length = length(unique(block)))
    ind <- vector(mode = "list", length = length(unique(block)))
    L <- vector(mode = "list", length = length(unique(block)))
    v_emp <- vector(mode = "list", length = length(unique(block)))
    
    p_emp <- v_vals_emp[,j]
    est_q0_pars <- fit.2g(q[which(p_emp>0.5)])$pars
    
    for(i in 1:length(unique(block))){
      k <- unique(block)[i]
      fold[[i]] <- which(block == k)
      ind[[i]] <- intersect(seq(1, length(p_emp), 1), fold[[i]])
      L[[i]] <- vl(p_emp, q, indices=ind[[i]], mode=2, fold=fold[[i]], gx = min(p_emp[ind[[i]]]));  # compute L curves
      v_emp[[i]] <- il(L[[i]], pi0_null = est_q0_pars[1], sigma_null = est_q0_pars[2], distx="norm") # integrate over L curves
      v_emp[[i]][which(v_emp[[i]]>1)] <- 1 # set any >1 equal to 1
    } 
    
    v_emp <- unlist(v_emp)
    v_vals_emp[,j+1] <- v_emp
    
    shared_CVs[[j]] <- shared_CV
    
  }
  
  # run BL
  BL_res <- lm_qvalue(p, X = q_vals[,c("q1","q2","q3","q4","q5")])
  
  list(data.frame(data, v_vals, v_vals_emp, q_vals, BL_qvals = BL_res$qvalues), shared_CVs)
}

# additional slurm commands
task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)

start <- 10*task_id + 1
end <- (10*task_id) + 5

out_res <- lapply(start:end, function(seed){
  set.seed(seed)
  func()
})

saveRDS(out_res, paste0("res/sims",task_id,".RDS"))
