# code to make bed file of asthma snps

x <- readRDS("asthma_snps.RDS")

chr <- paste0("chr",x$CHR19)
start <- x$BP19

bed <- data.frame(chr, start, start)

write.table(bed, "asthma_snps.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')