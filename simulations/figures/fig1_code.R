# code to make Figure 1

library(ggplot2)
library(cowplot)

#### collate simulation results

library(data.table)
library(dplyr)

files <- list.files(path = "../simA/res", pattern="*.RDS", full=TRUE)
data_indunif <- lapply(files, readRDS) # simulation A

files2 <- list.files(path = "../simB/res", pattern="*.RDS", full=TRUE)
data_relatedtrait <- lapply(files2, readRDS) # simulation B

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
# defined as those with r^2>=0.8 (r^2>=0.6) and r^2<=0.05 (r^2<=0.01) respectively
out_indunif <- lapply(data_indunif_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$assoc_0.6 <- as.numeric(x$max_r2>=0.6)
  x$noassoc_0.05 <- as.numeric(x$max_r2<=0.05)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

# sensitivity data.frame
sens_tmp <- lapply(out_indunif, function(x){
  
  sens_0.8 <- vector()
  sens_0.6 <- vector()
  sens_emp_0.8 <- vector()
  sens_emp_0.6 <- vector()
  
  for(i in 1:6){
    sens_0.8[i] <- length(which(x$assoc_0.8==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
    sens_0.6[i] <- length(which(x$assoc_0.6==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.6==1))
  }
  
  for(i in 1:6){
    sens_emp_0.8[i] <- length(which(x$assoc_0.8==1 & x[,12+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
    sens_emp_0.6[i] <- length(which(x$assoc_0.6==1 & x[,12+i]<=5*10^-6))/length(which(x$assoc_0.6==1))
  }
  
  list(sens_0.8, sens_0.6, sens_emp_0.8, sens_emp_0.6)
}) 

sens_0.8 <- lapply(sens_tmp, '[[', 1)
sens_0.8_func <- data.frame(do.call(rbind, sens_0.8))

sens_0.6 <- lapply(sens_tmp, '[[', 2)
sens_0.6_func <- data.frame(do.call(rbind, sens_0.6))

sens_emp_0.8 <- lapply(sens_tmp, '[[', 3)
sens_emp_0.8 <- data.frame(do.call(rbind, sens_emp_0.8))

sens_emp_0.6 <- lapply(sens_tmp, '[[', 4)
sens_emp_0.6 <- data.frame(do.call(rbind, sens_emp_0.6))

m_sens_0.8 <- melt(sens_0.8_func)
m_sens_0.6 <- melt(sens_0.6_func)
m_sens_emp_0.8 <- melt(sens_emp_0.8)
m_sens_emp_0.6 <- melt(sens_emp_0.6)

df <- data.frame(Iteration = c(m_sens_0.8$variable, m_sens_0.6$variable, m_sens_emp_0.8$variable, m_sens_emp_0.6$variable), Sensitivity = c(m_sens_0.8$value, m_sens_0.6$value, m_sens_emp_0.8$value, m_sens_emp_0.6$value), Method = c(rep("Flexible", times = length(c(m_sens_0.8$variable, m_sens_0.6$variable))), rep("Empirical", times = length(c(m_sens_emp_0.8$variable, m_sens_emp_0.6$variable)))), r2 = c(rep(0.8, length(m_sens_0.8$variable)), rep(0.6, length(m_sens_0.6$variable)), rep(0.8, length(m_sens_emp_0.8$variable)), rep(0.6, length(m_sens_emp_0.6$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Aleft <- ggplot(df, aes(x = Iteration, y = Sensitivity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8))+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+coord_cartesian(ylim = c(0.2,0.5)) + facet_grid(row = vars(r2)) + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

###########################

# specificity data.frame
specificity_tmp <- lapply(out_indunif, function(x){
  
  spec_0.05 <- vector()
  spec_0.01 <- vector()
  spec_emp_0.05 <- vector()
  spec_emp_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.05[i] <- length(which(x$noassoc_0.05==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.05==1))
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  for(i in 1:6){
    spec_emp_0.05[i] <- length(which(x$noassoc_0.05==1 & x[,12+i]>5*10^-6))/length(which(x$noassoc_0.05==1))
    spec_emp_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.05, spec_0.01, spec_emp_0.05, spec_emp_0.01)
}) 

spec_0.05 <- lapply(specificity_tmp, '[[', 1)
spec_0.05_func <- data.frame(do.call(rbind, spec_0.05))

spec_0.01 <- lapply(specificity_tmp, '[[', 2)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

spec_emp_0.05 <- lapply(specificity_tmp, '[[', 3)
spec_emp_0.05 <- data.frame(do.call(rbind, spec_emp_0.05))

spec_emp_0.01 <- lapply(specificity_tmp, '[[', 4)
spec_emp_0.01 <- data.frame(do.call(rbind, spec_emp_0.01))

m_spec_0.05 <- melt(spec_0.05_func)
m_spec_0.01 <- melt(spec_0.01_func)
m_spec_emp_0.05 <- melt(spec_emp_0.05)
m_spec_emp_0.01 <- melt(spec_emp_0.01)

df <- data.frame(Iteration = c(m_spec_0.05$variable, m_spec_0.01$variable, m_spec_emp_0.05$variable, m_spec_emp_0.01$variable), Specificity = c(m_spec_0.05$value, m_spec_0.01$value, m_spec_emp_0.05$value, m_spec_emp_0.01$value), Method = c(rep("Flexible", times = length(c(m_spec_0.05$variable, m_spec_0.01$variable))), rep("Empirical", times = length(c(m_spec_emp_0.05$variable, m_spec_emp_0.01$variable)))), r2 = c(rep(0.05, length(m_spec_0.05$variable)), rep(0.01, length(m_spec_0.01$variable)), rep(0.05, length(m_spec_emp_0.05$variable)), rep(0.01, length(m_spec_emp_0.01$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

df$r2_f = factor(df$r2, levels = c(0.05, 0.01))

Aright <- ggplot(df, aes(x = Iteration, y = Specificity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8))+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+ geom_hline(yintercept =  1-(5*10^-6), linetype = "dashed") + facet_grid(row = vars(r2_f)) + theme(panel.spacing = unit(1, "lines")) + coord_cartesian(ylim = c(0.9997, 1))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

top_panel <- plot_grid(Aleft + theme(legend.position = "none"), Aright + theme(legend.position = "none"), nrow = 1)

######################################################################
######################################################################
######################################################################

length(data_relatedtrait)
data_relatedtrait <- data_relatedtrait[1:100] # use 100 simulation results for figure

# convert p values to adjusted p values (using BH)
data_relatedtrait_qvals <- lapply(data_relatedtrait, function(x){
  
  for(i in 1:6){
    x[,6+i] <- p.adjust(x[,6+i], method = "BH")
  }
  
  for(i in 1:6){
    x[,12+i] <- p.adjust(x[,12+i], method = "BH")
  }
  
  x
})

out_relatedtrait <- lapply(data_relatedtrait_qvals, function(x){
  x$assoc_0.8 <- as.numeric(x$max_r2>=0.8)
  x$assoc_0.6 <- as.numeric(x$max_r2>=0.6)
  x$noassoc_0.05 <- as.numeric(x$max_r2<=0.05)
  x$noassoc_0.01 <- as.numeric(x$max_r2<=0.01)
  x
})

sens_tmp <- lapply(out_relatedtrait, function(x){
  
  sens_0.8 <- vector()
  sens_0.6 <- vector()
  sens_emp_0.8 <- vector()
  sens_emp_0.6 <- vector()
  
  for(i in 1:6){
    sens_0.8[i] <- length(which(x$assoc_0.8==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
    sens_0.6[i] <- length(which(x$assoc_0.6==1 & x[,6+i]<=5*10^-6))/length(which(x$assoc_0.6==1))
  }
  
  for(i in 1:6){
    sens_emp_0.8[i] <- length(which(x$assoc_0.8==1 & x[,12+i]<=5*10^-6))/length(which(x$assoc_0.8==1))
    sens_emp_0.6[i] <- length(which(x$assoc_0.6==1 & x[,12+i]<=5*10^-6))/length(which(x$assoc_0.6==1))
  }
  
  list(sens_0.8, sens_0.6, sens_emp_0.8, sens_emp_0.6)
}) 

sens_0.8 <- lapply(sens_tmp, '[[', 1)
sens_0.8_func <- data.frame(do.call(rbind, sens_0.8))

sens_0.6 <- lapply(sens_tmp, '[[', 2)
sens_0.6_func <- data.frame(do.call(rbind, sens_0.6))

sens_emp_0.8 <- lapply(sens_tmp, '[[', 3)
sens_emp_0.8 <- data.frame(do.call(rbind, sens_emp_0.8))

sens_emp_0.6 <- lapply(sens_tmp, '[[', 4)
sens_emp_0.6 <- data.frame(do.call(rbind, sens_emp_0.6))

m_sens_0.8 <- melt(sens_0.8_func)
m_sens_0.6 <- melt(sens_0.6_func)
m_sens_emp_0.8 <- melt(sens_emp_0.8)
m_sens_emp_0.6 <- melt(sens_emp_0.6)

df <- data.frame(Iteration = c(m_sens_0.8$variable, m_sens_0.6$variable, m_sens_emp_0.8$variable, m_sens_emp_0.6$variable), Sensitivity = c(m_sens_0.8$value, m_sens_0.6$value, m_sens_emp_0.8$value, m_sens_emp_0.6$value), Method = c(rep("Flexible", times = length(c(m_sens_0.8$variable, m_sens_0.6$variable))), rep("Empirical", times = length(c(m_sens_emp_0.8$variable, m_sens_emp_0.6$variable)))), r2 = c(rep(0.8, length(m_sens_0.8$variable)), rep(0.6, length(m_sens_0.6$variable)), rep(0.8, length(m_sens_emp_0.8$variable)), rep(0.6, length(m_sens_emp_0.6$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

Bleft <- ggplot(df, aes(x = Iteration, y = Sensitivity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8))+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+coord_cartesian(ylim = c(0.2,0.5)) + facet_grid(row = vars(r2)) + theme(panel.spacing = unit(1, "lines"))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

###########################

specificity_tmp <- lapply(out_relatedtrait, function(x){
  
  spec_0.05 <- vector()
  spec_0.01 <- vector()
  spec_emp_0.05 <- vector()
  spec_emp_0.01 <- vector()
  
  for(i in 1:6){
    spec_0.05[i] <- length(which(x$noassoc_0.05==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.05==1))
    spec_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,6+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  for(i in 1:6){
    spec_emp_0.05[i] <- length(which(x$noassoc_0.05==1 & x[,12+i]>5*10^-6))/length(which(x$noassoc_0.05==1))
    spec_emp_0.01[i] <- length(which(x$noassoc_0.01==1 & x[,12+i]>5*10^-6))/length(which(x$noassoc_0.01==1))
  }
  
  list(spec_0.05, spec_0.01, spec_emp_0.05, spec_emp_0.01)
}) 

spec_0.05 <- lapply(specificity_tmp, '[[', 1)
spec_0.05_func <- data.frame(do.call(rbind, spec_0.05))

spec_0.01 <- lapply(specificity_tmp, '[[', 2)
spec_0.01_func <- data.frame(do.call(rbind, spec_0.01))

spec_emp_0.05 <- lapply(specificity_tmp, '[[', 3)
spec_emp_0.05 <- data.frame(do.call(rbind, spec_emp_0.05))

spec_emp_0.01 <- lapply(specificity_tmp, '[[', 4)
spec_emp_0.01 <- data.frame(do.call(rbind, spec_emp_0.01))

m_spec_0.05 <- melt(spec_0.05_func)
m_spec_0.01 <- melt(spec_0.01_func)
m_spec_emp_0.05 <- melt(spec_emp_0.05)
m_spec_emp_0.01 <- melt(spec_emp_0.01)

df <- data.frame(Iteration = c(m_spec_0.05$variable, m_spec_0.01$variable, m_spec_emp_0.05$variable, m_spec_emp_0.01$variable), Specificity = c(m_spec_0.05$value, m_spec_0.01$value, m_spec_emp_0.05$value, m_spec_emp_0.01$value), Method = c(rep("Flexible", times = length(c(m_spec_0.05$variable, m_spec_0.01$variable))), rep("Empirical", times = length(c(m_spec_emp_0.05$variable, m_spec_emp_0.01$variable)))), r2 = c(rep(0.05, length(m_spec_0.05$variable)), rep(0.01, length(m_spec_0.01$variable)), rep(0.05, length(m_spec_emp_0.05$variable)), rep(0.01, length(m_spec_emp_0.01$variable))))

df$Iteration[which(df$Iteration==1)] <- 0
df$Iteration[which(df$Iteration==2)] <- 1
df$Iteration[which(df$Iteration==3)] <- 2
df$Iteration[which(df$Iteration==4)] <- 3
df$Iteration[which(df$Iteration==5)] <- 4
df$Iteration[which(df$Iteration==6)] <- 5

df$Iteration <- as.factor(df$Iteration)

df$r2_f = factor(df$r2, levels = c(0.05, 0.01))

Bright <- ggplot(df, aes(x = Iteration, y = Specificity, col = Method)) + geom_point(stat = "summary", fun = "mean", position = position_dodge(0.8))+ 
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1), position = position_dodge(0.8), width = 0.6) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + scale_colour_manual(values = c("orchid4","steelblue1"))+ geom_hline(yintercept =  1-(5*10^-6), linetype = "dashed") + facet_grid(row = vars(r2_f)) + theme(panel.spacing = unit(1, "lines"))+ coord_cartesian(ylim = c(0.9997, 1))+ theme(text = element_text(size = 10), axis.title = element_text(size = 11))

bottom_panel <- plot_grid(Bleft + theme(legend.position = "none"), Bright + theme(legend.position = "none"), nrow = 1)

final <- plot_grid(top_panel, bottom_panel, nrow = 2, labels = c("A","B"))

legend <- get_legend(
  # create some space to the left of the legend
  Bleft + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p = plot_grid(final, legend, rel_widths = c(3, .6))

ggsave("fig1.png", p, dev = "png", width = 15, height = 21, units = "cm")