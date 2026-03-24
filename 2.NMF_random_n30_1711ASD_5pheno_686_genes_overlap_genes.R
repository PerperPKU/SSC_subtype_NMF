## This script is to choose a optimal rank for NMF analysis.
## We performed 30 runs to determine the optimal factorization rank
## Installing the foreach and doParallel packages alongside NMF is strongly recommended, 
## as they enable transparent parallel computation across all available CPU cores, 
## significantly accelerating the analysis.
## plot(nmf.results) will produce supplementary figure 12 
library(NMF)
expr=read.table("dir/results/1711_ASDs_gene_counts_normed_by_limma_686genes_matrix.txt",header=T)
print(dim(expr))
expr<-data.frame(expr)
rank<-c(2:10)
nmf.results<- nmf(expr,rank,nrun=30)
save(nmf.results,file="dir/results/1711_686genes_NMF_results_rank_2_10_nrun50_bootstraped_overlapgenes_asd.Rdata")
pdf('dir/results/1711_686genes_NMF_results_rank_2_10_nrun50_bootstraped_overlapgenes_asd.pdf')
plot(nmf.results)
dev.off()

