##### code for dependent continous q (simulation D)

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
  q_vals <- data.frame("q1" = runif(length(p), 0, 1), "q2" = runif(length(p), 0, 1), "q3" = runif(length(p), 0, 1), "q4" = runif(length(p), 0, 1), "q5" = runif(length(p), 0, 1))
  
  # ensure random seeds
  s <- sample(1:100000, 1)
  niter <- 5
  
  for(j in 1:niter){
    
    set.seed(j*s)
    
    p <- v_vals[,j]
    n <- length(p)
    
    # simulate functional and non-functionals from seperate distributions
    z <- runif(n)
    weight <- sample(c(0.6,0.7,0.8,0.9,0.95), 1)
    mean1 <- sample(c(2.5, 3, 4), 1)
    mean2 <- sample(c(-3, -2, -1.5), 1)
    
    mixture_comp1 <- function(x) rnorm(x, mean = mean1, sd = 1)
    mixture_comp2 <- function(x) rnorm(x, mean = mean2, sd = 0.5)
    
    q <- ifelse(z<weight, mixture_comp1(n), mixture_comp2(n))
    
    q[which(functionals==1)] <- ifelse(z[which(functionals==1)]<weight, mixture_comp2(length(which(functionals==1))), mixture_comp1(length(which(functionals==1))))
    
    # ensure that the sign of the correlation is the same for 
    # the whole data set and the independent subset 
    while( sign(cor(p[indep_index], q[indep_index], method="spearman"))!= sign(cor(p, q, method="spearman")) ){
      z <- runif(n)
      weight <- sample(c(0.6,0.7,0.8,0.9,0.95), 1)
      
      # sample non-functionals from 0.9*MC1 + 0.1*MC2
      q <- ifelse(z<weight, mixture_comp1(n), mixture_comp2(n))
      
      q[which(functionals==1)] <- ifelse(z[which(functionals==1)]<weight, mixture_comp2(length(which(functionals==1))), mixture_comp1(length(which(functionals==1))))
    } 
    
    q_vals[,j] <- q
    
    print(paste0("Iteration: ", j))
    print(paste0("Correlation is: ",cor(p, q)))
    print(paste0("Ind correlation is: ",cor(p[indep_index], q[indep_index])))
    print(paste0("q is :", c(q[1], q[2])))
    print(paste0("weight is: ", weight))
    print(paste0("means are: ", c(mean1, mean2)))
    
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