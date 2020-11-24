source("fig_source.R")

# code to make Figure 3

library(ggplot2)
library(cowplot)

#### collate simulation results

library(data.table)
library(dplyr)

# simulation E
files <- list.files(path = "../simE/res", pattern="*.RDS", full=TRUE)
data_sameq <- lapply(files, readRDS)

length(data_sameq)
data_sameq <- data_sameq[1:100] # use 100 simulation results for figure

data_sameq_qvals <- lapply(data_sameq, function(x){
  
  for(i in 1:6){
    x[,6+i] <- p.adjust(x[,6+i], method = "BH")
  }
  
  x
})

out_sameq <- lapply(data_sameq_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

sens_tmp <- lapply(out_sameq, function(x){
  
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

left <- ggplot(df, aes(x = Iteration, y = Sensitivity)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8))+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") +coord_cartesian(ylim = c(0.25,0.8))  + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

#############

specificity_tmp <- lapply(out_sameq, function(x){
  
  spec_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.01)
}) 

spec_0.01 <- lapply(specificity_tmp, '[[', 1)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

m_spec_0.01 <- melt(spec_0.01_func)

df <- data.frame(Iteration = c(m_spec_0.01$variable), Specificity = c( m_spec_0.01$value), r2 = c(rep(0.01, times = length(m_spec_0.01$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

right <- ggplot(df, aes(x = Iteration, y = Specificity)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8))+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.4) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_hline(yintercept =  1-(5*10^-6), linetype = "dashed")+ theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

p = plot_grid(left + theme(legend.position = "none"), right, nrow = 1, labels = c("A","B"))

ggsave("fig2.png", plot = p, dev = "png", width = 15, height = 10, units = "cm")

