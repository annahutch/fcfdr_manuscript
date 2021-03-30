# code to make Figure 1

library(ggplot2)
library(cowplot)

#### collate simulation results

library(data.table)
library(dplyr)

# simulation A
files <- list.files(path = "../simA/res", pattern="*.RDS", full=TRUE)
data_indunif <- lapply(files, readRDS) 

# simulation B
files2 <- list.files(path = "../simB/res", pattern="*.RDS", full=TRUE)
data_relatedtrait <- lapply(files2, readRDS) 

# simulation C
files3 <- list.files(path = "../simC/res", pattern="*.RDS", full=TRUE)
data_indepcont <- lapply(files3, readRDS) 

# simulation D
files4 <- list.files(path = "../simD/res", pattern="*.RDS", full=TRUE)
data_probfunc <- lapply(files4, readRDS) 

#################

no_CV <- lapply(data_indunif, function(x) length(which(x$CV==1))) %>% unlist()

length(data_indunif)
data_indunif <- data_indunif[1:100] # use 100 simulation results for figure

# convert p values to adjusted p values (using BH)
data_indunif_qvals <- lapply(data_indunif, function(x){
  
  for(i in 1:6){
    x[,6+i] <- p.adjust(x[,6+i], method = "BH")
  }
  
  for(i in 1:6){
    x[,12+i] <- p.adjust(x[,12+i], method = "BH")
  }
  
  x
})

# find associated and non-associated sets of SNPs
# defined as those with r^2>=0.8 and r^2<=0.01 respectively
out_indunif <- lapply(data_indunif_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

# sensitivity data.frame
sens_tmp <- lapply(out_indunif, function(x){
  
  sens_0.8 <- vector()
  sens_emp_0.8 <- vector()
  
  for(i in 1:6){
    sens_0.8[i] <- length(which(x$assoc_0.8==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
  }
  
  for(i in 1:6){
    sens_emp_0.8[i] <- length(which(x$assoc_0.8==1 & x[,12+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
  }
  
  list(sens_0.8, sens_emp_0.8)
}) 

sens_0.8 <- lapply(sens_tmp, '[[', 1)
sens_0.8_func <- data.frame(do.call(rbind, sens_0.8))

sens_emp_0.8 <- lapply(sens_tmp, '[[', 2)
sens_emp_0.8 <- data.frame(do.call(rbind, sens_emp_0.8))

m_sens_0.8 <- melt(sens_0.8_func)
m_sens_emp_0.8 <- melt(sens_emp_0.8)

df <- data.frame(Iteration = c(m_sens_0.8$variable, m_sens_emp_0.8$variable), Sensitivity = c(m_sens_0.8$value, m_sens_emp_0.8$value), Method = c(rep("Flexible", times = length(c(m_sens_0.8$variable))), rep("Empirical", times = length(c(m_sens_emp_0.8$variable)))), r2 = c(rep(0.8, length(c(m_sens_0.8$variable, m_sens_emp_0.8$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Aleft <- ggplot(df, aes(x = Iteration, y = Sensitivity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+coord_cartesian(ylim = c(0.25,0.5)) + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

###########################

# specificity data.frame
specificity_tmp <- lapply(out_indunif, function(x){
  
  spec_0.01 <- vector()
  spec_emp_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  for(i in 1:6){
    spec_emp_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.01,spec_emp_0.01)
}) 


spec_0.01 <- lapply(specificity_tmp, '[[', 1)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

spec_emp_0.01 <- lapply(specificity_tmp, '[[', 2)
spec_emp_0.01 <- data.frame(do.call(rbind, spec_emp_0.01))

m_spec_0.01 <- melt(spec_0.01_func)
m_spec_emp_0.01 <- melt(spec_emp_0.01)

df <- data.frame(Iteration = c(m_spec_0.01$variable, m_spec_emp_0.01$variable), Specificity = c(m_spec_0.01$value, m_spec_emp_0.01$value), Method = c(rep("Flexible", times = length(c(m_spec_0.01$variable))), rep("Empirical", times = length(c(m_spec_emp_0.01$variable)))), r2 = c(rep(0.01, length(c(m_spec_0.01$variable, m_spec_emp_0.01$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Amid <-  ggplot(df, aes(x = Iteration, y = Specificity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), size = 0.8) + geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1")) + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0.9999, 1))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ scale_y_continuous(breaks = c(0.9999, 1))

###########################

# fdr data.frame
fdr_tmp <- lapply(out_indunif, function(x){
  
  fdr_func <- vector()
  fdr_emp <- vector()
  
  for(i in 1:6){
    fdr_func[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]<=0.05))/length(which(x[,6+i]<=0.05))
  }
  
  for(i in 1:6){
    fdr_emp[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]<=0.05))/length(which(x[,12+i]<=0.05))
  }
  
  list(fdr_func, fdr_emp)
})

fdr_func_tmp <- lapply(fdr_tmp, '[[', 1)
fdr_func_final <- data.frame(do.call(rbind, fdr_func_tmp))

fdr_emp_tmp <- lapply(fdr_tmp, '[[', 2)
fdr_emp_final <- data.frame(do.call(rbind, fdr_emp_tmp))

fdr_func_m <- melt(fdr_func_final)
fdr_emp_m <- melt(fdr_emp_final)

df <- data.frame(Iteration = c(fdr_func_m$variable, fdr_emp_m$variable), FDR = c(fdr_func_m$value, fdr_emp_m$value), Method = c(rep("Flexible", times = length(c(fdr_func_m$variable))), rep("Empirical", times = length(c(fdr_emp_m$variable)))), r2 = c(rep(0.01, length(c(fdr_func_m$variable, fdr_emp_m$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Aright <- ggplot(df, aes(x = Iteration, y = FDR, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+ geom_hline(yintercept =  0.05, linetype = "dashed") + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0, 0.05))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ scale_y_continuous(breaks = c(0, 0.05))

A_final <- plot_grid(Aleft + theme(legend.position = "none"), Amid + theme(legend.position = "none"), Aright + theme(legend.position = "none"), nrow = 1)

######################################################################
######################################################################
######################################################################

length(data_relatedtrait)
data_relatedtrait <- data_relatedtrait[1:100] # use 100 simulation results for figure

# convert p values to adjusted p values (using BH)
data_related_trait_qvals <- lapply(data_relatedtrait, function(x){
  
  for(i in 1:6){
    x[,6+i] <- p.adjust(x[,6+i], method = "BH")
  }
  
  for(i in 1:6){
    x[,12+i] <- p.adjust(x[,12+i], method = "BH")
  }
  
  x
})

out_related_trait <- lapply(data_related_trait_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

sens_tmp <- lapply(out_related_trait, function(x){
  
  sens_0.8 <- vector()
  sens_emp_0.8 <- vector()
  
  for(i in 1:6){
    sens_0.8[i] <- length(which(x$assoc_0.8==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
  }
  
  for(i in 1:6){
    sens_emp_0.8[i] <- length(which(x$assoc_0.8==1 & x[,12+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
  }
  
  list(sens_0.8, sens_emp_0.8)
}) 

sens_0.8 <- lapply(sens_tmp, '[[', 1)
sens_0.8_func <- data.frame(do.call(rbind, sens_0.8))

sens_emp_0.8 <- lapply(sens_tmp, '[[', 2)
sens_emp_0.8 <- data.frame(do.call(rbind, sens_emp_0.8))

m_sens_0.8 <- melt(sens_0.8_func)
m_sens_emp_0.8 <- melt(sens_emp_0.8)

df <- data.frame(Iteration = c(m_sens_0.8$variable, m_sens_emp_0.8$variable), Sensitivity = c(m_sens_0.8$value, m_sens_emp_0.8$value), Method = c(rep("Flexible", times = length(c(m_sens_0.8$variable))), rep("Empirical", times = length(c(m_sens_emp_0.8$variable)))), r2 = c(rep(0.8, length(c(m_sens_0.8$variable, m_sens_emp_0.8$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Bleft <- ggplot(df, aes(x = Iteration, y = Sensitivity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+coord_cartesian(ylim = c(0.25,0.5)) + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

###########################

specificity_tmp <- lapply(out_related_trait, function(x){
  
  spec_0.01 <- vector()
  spec_emp_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  for(i in 1:6){
    spec_emp_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.01,spec_emp_0.01)
}) 


spec_0.01 <- lapply(specificity_tmp, '[[', 1)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

spec_emp_0.01 <- lapply(specificity_tmp, '[[', 2)
spec_emp_0.01 <- data.frame(do.call(rbind, spec_emp_0.01))

m_spec_0.01 <- melt(spec_0.01_func)
m_spec_emp_0.01 <- melt(spec_emp_0.01)

df <- data.frame(Iteration = c(m_spec_0.01$variable, m_spec_emp_0.01$variable), Specificity = c(m_spec_0.01$value, m_spec_emp_0.01$value), Method = c(rep("Flexible", times = length(c(m_spec_0.01$variable))), rep("Empirical", times = length(c(m_spec_emp_0.01$variable)))), r2 = c(rep(0.01, length(c(m_spec_0.01$variable, m_spec_emp_0.01$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Bmid <- ggplot(df, aes(x = Iteration, y = Specificity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), size = 0.8) + geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1")) + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0.9999, 1))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ scale_y_continuous(breaks = c(0.9999, 1))

###########################

fdr_tmp <- lapply(out_related_trait, function(x){
  
  fdr_func <- vector()
  fdr_emp <- vector()
  
  for(i in 1:6){
    fdr_func[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]<=0.05))/length(which(x[,6+i]<=0.05))
  }
  
  for(i in 1:6){
    fdr_emp[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]<=0.05))/length(which(x[,12+i]<=0.05))
  }
  
  list(fdr_func, fdr_emp)
})

fdr_func_tmp <- lapply(fdr_tmp, '[[', 1)
fdr_func_final <- data.frame(do.call(rbind, fdr_func_tmp))

fdr_emp_tmp <- lapply(fdr_tmp, '[[', 2)
fdr_emp_final <- data.frame(do.call(rbind, fdr_emp_tmp))

fdr_func_m <- melt(fdr_func_final)
fdr_emp_m <- melt(fdr_emp_final)

df <- data.frame(Iteration = c(fdr_func_m$variable, fdr_emp_m$variable), FDR = c(fdr_func_m$value, fdr_emp_m$value), Method = c(rep("Flexible", times = length(c(fdr_func_m$variable))), rep("Empirical", times = length(c(fdr_emp_m$variable)))), r2 = c(rep(0.01, length(c(fdr_func_m$variable, fdr_emp_m$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Bright <- ggplot(df, aes(x = Iteration, y = FDR, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+ geom_hline(yintercept =  0.05, linetype = "dashed") + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0, 0.2))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ scale_y_continuous(breaks = c(0, 0.1, 0.2))

B_final <- plot_grid(Bleft + theme(legend.position = "none"), Bmid + theme(legend.position = "none"), Bright + theme(legend.position = "none"), nrow = 1)

###########################

AB <- plot_grid(A_final, B_final, nrow = 2, labels = c("A","B"))

legend <- get_legend(
  # create some space to the left of the legend
  Bleft + theme(legend.box.margin = margin(0, 0, 0, 12))
)

AB_final = plot_grid(AB, legend, rel_widths = c(3, .6))

###############################

length(data_indepcont)
data_indepcont <- data_indepcont[1:100]

data_indepcont_qvals <- lapply(data_indepcont, function(x){
  
  for(i in 1:6){
    x[,6+i] <- p.adjust(x[,6+i], method = "BH")
  }
  
  x
})

out_indepcont <- lapply(data_indepcont_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

sens_tmp <- lapply(out_indepcont, function(x){
  
  sens_0.8 <- vector()
  
  for(i in 1:6){
    sens_0.8[i] <- length(which(x$assoc_0.8==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
  }
  
  list(sens_0.8)
}) 

sens_0.8 <- lapply(sens_tmp, '[[', 1)
sens_0.8_func <- data.frame(do.call(rbind, sens_0.8))

m_sens_0.8 <- melt(sens_0.8_func)

df <- data.frame(Iteration = c(m_sens_0.8$variable), Sensitivity = m_sens_0.8$value, r2 = c(rep(0.8, times = length(m_sens_0.8$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Cleft <- ggplot(df, aes(x = Iteration, y = Sensitivity)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), col = "steelblue1", size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6, col = "steelblue1") + theme_cowplot(12) + background_grid(major = "xy", minor = "none") +coord_cartesian(ylim = c(0.25,0.5)) + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

#############

specificity_tmp <- lapply(out_indepcont, function(x){
  
  spec_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.01)
}) 

spec_0.01 <- lapply(specificity_tmp, '[[', 1)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

m_spec_0.01 <- melt(spec_0.01_func)

df <- data.frame(Iteration = c(m_spec_0.01$variable), Specificity = m_spec_0.01$value, r2 = c(rep(0.01, times = length(m_spec_0.01$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Cmid <- ggplot(df, aes(x = Iteration, y = Specificity)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), col = "steelblue1", size = 0.8) + geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.4, col = "steelblue1") + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ coord_cartesian(ylim = c(0.9999, 1))+ scale_y_continuous(breaks = c(0.9999, 1))

#############

fdr_tmp <- lapply(out_indepcont, function(x){
  
  fdr_func <- vector()
  fdr_emp <- vector()
  
  for(i in 1:6){
    fdr_func[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]<=0.05))/length(which(x[,6+i]<=0.05))
  }
  
  for(i in 1:6){
    fdr_emp[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]<=0.05))/length(which(x[,12+i]<=0.05))
  }
  
  list(fdr_func, fdr_emp)
})

fdr_func_tmp <- lapply(fdr_tmp, '[[', 1)
fdr_func_final <- data.frame(do.call(rbind, fdr_func_tmp))

fdr_emp_tmp <- lapply(fdr_tmp, '[[', 2)
fdr_emp_final <- data.frame(do.call(rbind, fdr_emp_tmp))

fdr_func_m <- melt(fdr_func_final)
fdr_emp_m <- melt(fdr_emp_final)

df <- data.frame(Iteration = c(fdr_func_m$variable, fdr_emp_m$variable), FDR = c(fdr_func_m$value, fdr_emp_m$value), Method = c(rep("Flexible", times = length(c(fdr_func_m$variable))), rep("Empirical", times = length(c(fdr_emp_m$variable)))), r2 = c(rep(0.01, length(c(fdr_func_m$variable, fdr_emp_m$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Cright <- ggplot(df, aes(x = Iteration, y = FDR, col = Method)) + geom_point(stat = "summary", fun = "mean", size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+ geom_hline(yintercept =  0.05, linetype = "dashed") + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0, 0.05))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ scale_y_continuous(breaks = c(0, 0.05))

C_final <- plot_grid(Cleft + theme(legend.position = "none"), Cmid + theme(legend.position = "none"), Cright + theme(legend.position = "none"), nrow = 1)

######################################################################
######################################################################
######################################################################

length(data_probfunc)
data_probfunc <- data_probfunc[1:100]

data_probfunc_qvals <- lapply(data_probfunc, function(x){
  
  for(i in 1:6){
    x[,6+i] <- p.adjust(x[,6+i], method = "BH")
  }
  
  x
})

out_probfunc <- lapply(data_probfunc_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

sens_tmp <- lapply(out_probfunc, function(x){
  
  sens_0.8 <- vector()
  
  for(i in 1:6){
    sens_0.8[i] <- length(which(x$assoc_0.8==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
  }
  
  list(sens_0.8)
}) 

sens_0.8 <- lapply(sens_tmp, '[[', 1)
sens_0.8_func <- data.frame(do.call(rbind, sens_0.8))

m_sens_0.8 <- melt(sens_0.8_func)

df <- data.frame(Iteration = c(m_sens_0.8$variable), Sensitivity = c(m_sens_0.8$value), r2 = c(rep(0.8, times = length(m_sens_0.8$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Dleft <- ggplot(df, aes(x = Iteration, y = Sensitivity)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), col = "steelblue1", size = 0.8)+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6, col = "steelblue1") + theme_cowplot(12) + background_grid(major = "xy", minor = "none") +coord_cartesian(ylim = c(0.25,0.5))+ theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

#############

specificity_tmp <- lapply(out_probfunc, function(x){
  
  spec_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.01)
}) 

spec_0.01 <- lapply(specificity_tmp, '[[', 1)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

m_spec_0.01 <- melt(spec_0.01_func)

df <- data.frame(Iteration = c(m_spec_0.01$variable), Specificity = c(m_spec_0.01$value), r2 = c(rep(0.01, times = length(m_spec_0.01$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Dmid <- ggplot(df, aes(x = Iteration, y = Specificity)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8), col = "steelblue1", size = 0.8) + geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.4, col = "steelblue1") + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ coord_cartesian(ylim = c(0.9999, 1))+ scale_y_continuous(breaks = c(0.9999, 1))

#############

fdr_tmp <- lapply(out_probfunc, function(x){
  
  fdr_func <- vector()
  fdr_emp <- vector()
  
  for(i in 1:6){
    fdr_func[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]<=0.05))/length(which(x[,6+i]<=0.05))
  }
  
  for(i in 1:6){
    fdr_emp[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]<=0.05))/length(which(x[,12+i]<=0.05))
  }
  
  list(fdr_func, fdr_emp)
})

fdr_func_tmp <- lapply(fdr_tmp, '[[', 1)
fdr_func_final <- data.frame(do.call(rbind, fdr_func_tmp))

fdr_emp_tmp <- lapply(fdr_tmp, '[[', 2)
fdr_emp_final <- data.frame(do.call(rbind, fdr_emp_tmp))

fdr_func_m <- melt(fdr_func_final)
fdr_emp_m <- melt(fdr_emp_final)

df <- data.frame(Iteration = c(fdr_func_m$variable, fdr_emp_m$variable), FDR = c(fdr_func_m$value, fdr_emp_m$value), Method = c(rep("Flexible", times = length(c(fdr_func_m$variable))), rep("Empirical", times = length(c(fdr_emp_m$variable)))), r2 = c(rep(0.01, length(c(fdr_func_m$variable, fdr_emp_m$variable)))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Dright <- ggplot(df, aes(x = Iteration, y = FDR, col = Method)) + geom_point(stat = "summary", fun = "mean", size = 0.8)+
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+ geom_hline(yintercept =  0.05, linetype = "dashed") + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0, 0.05))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))+ scale_y_continuous(breaks = c(0, 0.05))

D_final <- plot_grid(Dleft + theme(legend.position = "none"), Dmid + theme(legend.position = "none"), Dright + theme(legend.position = "none"), nrow = 1)

CD <- plot_grid(C_final, D_final, nrow = 2, labels = c("C","D"))

legend <- get_legend(
  # create some space to the left of the legend
  Bleft + theme(legend.position = "none") + theme(legend.box.margin = margin(0, 0, 0, 12))
)

CD_final = plot_grid(CD, legend, rel_widths = c(3, .6))

plot_final <- plot_grid(AB_final, CD_final, nrow = 2)

ggsave("fig1.png", plot = plot_final, dev = "png", units = "cm", width = 15, height = 21)
