# code to collate H3K27ac fold change values for asthma SNPs into a single .RDS file

library(data.table)
library(dplyr)

# randomly select one to keep
# https://amywhiteheadresearch.wordpress.com/2013/01/22/randomly-deleting-duplicate-rows-from-a-dataframe-2/
duplicated.random = function(x, incomparables = FALSE, ...) 
{ 
  if ( is.vector(x) ) 
  { 
    permutation = sample(length(x)) 
    x.perm      = x[permutation] 
    result.perm = duplicated(x.perm, incomparables, ...) 
    result      = result.perm[order(permutation)] 
    return(result) 
  } 
  else if ( is.matrix(x) ) 
  { 
    permutation = sample(nrow(x)) 
    x.perm      = x[permutation,] 
    result.perm = duplicated(x.perm, incomparables, ...) 
    result      = result.perm[order(permutation)] 
    return(result) 
  } 
  else 
  { 
    stop(paste("duplicated.random() only supports vectors", 
               "matrices for now.")) 
  } 
} 

files <- list.files(path = "asthmasnp_match", pattern="", full=TRUE)

# randomly select rows where a SNP matches twice (i.e. falls on a boundary)
# to obtain the H3K27ac fold change values for each SNP in each cell type
res <- lapply(files, function(z){
  
  y <- read.table(z)
  print(dim(y)[1]>1968651)
  
  dups <- duplicated.random(as.matrix(y[,c(1:3)]))
  
  # final data set with randomly kept counts for SNPs that match twice
  final <- y[which(!dups),c("V8")]
  
  final
})

df <- t(data.frame(matrix(unlist(res), nrow=length(res), byrow=T)))

names <- sapply(files, function(x) substr(x, 17, 20)) %>% unname()
colnames(df) <- names
rownames(df) <- NULL

## rename columns with more intuitive names

colnames(df) <- gsub("E029", "CD14_Primary_Cells", colnames(df))
colnames(df) <- gsub("E045", "CD4+_CD25int_CD127+_Tmem_Primary_Cells", colnames(df))
colnames(df) <- gsub("E044", "CD4+_CD25+_CD127-_Treg_Primary_Cells", colnames(df))
colnames(df) <- gsub("E043", "CD4+_CD25-_Th_Primary_Cells", colnames(df))
colnames(df) <- gsub("E039", "CD4+_CD25-_CD45RA+_Naive_Primary_Cells", colnames(df))
colnames(df) <- gsub("E040", "CD4+_CD25-_CD45RO+_Memory_Primary_Cells", colnames(df))
colnames(df) <- gsub("E037", "CD4_Memory_Primary_Cells", colnames(df))
colnames(df) <- gsub("E048", "CD8_Memory_Primary_Cells", colnames(df))
colnames(df) <- gsub("E038", "CD4_Naive_Primary_Cells", colnames(df))
colnames(df) <- gsub("E047", "CD8_Naive_Primary_Cells", colnames(df))
colnames(df) <- gsub("E032", "CD19_Primary_Cells_Peripheral_UW", colnames(df))
colnames(df) <- gsub("E096", "Lung", colnames(df))
colnames(df) <- gsub("E046", "CD56_Primary_Cells", colnames(df))

# take logs
df <- log(df+1)

# only use H3K27ac
q1 <- rowMeans(df[,c(2,3,4,5,6,7,8,9,10,11,12)])
q2 <- rowMeans(df[,c(1,13)])

saveRDS(data.frame(q1, q2), "H3K27ac_qs.RDS")