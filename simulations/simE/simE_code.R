##### code for repeatedly simulating over q sampled from the same distribution (simulation E)

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
  
  # to make dependent auxiliary data we need to define functional SNPs
  # define these as CVs +/- 10000 bp
  CV_pos <- data$pos[which(data$CV==1)]
  intervals <- sapply(CV_pos, function(x) c(x-10000,x+10000), simplify=FALSE) %>% unlist()
  functionals <- as.integer(inrange(data$pos, lower = intervals[c(TRUE,FALSE)], upper = intervals[c(FALSE,TRUE)]))
  
  # make containers for V values (with first column original P values for first iteration)
  v_vals <- data.frame("p" = p, "v1" = p, "v2" = p, "v3" = p, "v4" = p, "v5" = p)
  
  # ensure random seeds
  s <- sample(1:100000, 1)
  niter <- 5
  
  # use q sampled from the same distribution
  n <- length(p)
  q1 <- rnorm(n, mean = 3, sd = 2)
  q1[which(functionals==1)] <- rnorm(length(which(functionals==1)), mean = -2, sd = 0.5)
  
  q2 <- rnorm(n, mean = 3, sd = 2)
  q2[which(functionals==1)] <- rnorm(length(which(functionals==1)), mean = -2, sd = 0.5)
  
  q3 <- rnorm(n, mean = 3, sd = 2)
  q3[which(functionals==1)] <- rnorm(length(which(functionals==1)), mean = -2, sd = 0.5)
  
  q4 <- rnorm(n, mean = 3, sd = 2)
  q4[which(functionals==1)] <- rnorm(length(which(functionals==1)), mean = -2, sd = 0.5)
  
  q5 <- rnorm(n, mean = 3, sd = 2)
  q5[which(functionals==1)] <- rnorm(length(which(functionals==1)), mean = -2, sd = 0.5)
  
  q_vals <- data.frame("q1" = q1, "q2" = q2, "q3" = q3, "q4" = q4, "q5" = q5)
  
  for(j in 1:niter){
    
    set.seed(j*s)
    
    p <- v_vals[,j]
    
    q <- q_vals[,j]
    
    print(paste0("Iteration: ", j))
    print(paste0("Correlation is: ",cor(p, q)))
    print(paste0("Ind correlation is: ",cor(p[indep_index], q[indep_index])))
    print(paste0("q is :", c(q[1], q[2])))
    print(paste0("seed is: ", j*s))
    
    # run flexible cFDR
    res_full <- flexible_cfdr(p = p, q = q, indep_index, gridp = 20)
    
    res <- res_full[[1]]
    v_vals[,j+1] <- res$v
    
  }
  
  data.frame(data, v_vals, q_vals)
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