## This script is to perform NMF analysis, 
## and nmf.results include w matrix and h matirx, helping us cluster samples and analyze metagenes.
start_time <- Sys.time()
library(NMF)
library(data.table)
expr=read.table("dir/results/1711_ASDs_gene_counts_normed_by_limma_686genes_matrix.txt",header=T)
print(dim(expr))
expr<-data.frame(expr)
rank=3
nmf.results<- nmf(expr,rank,nrun=300)
nmf.fit<-fit(nmf.results)
w <- basis(nmf.results)
h <- coef(nmf.results)
feature_score <- featureScore(nmf.results)
save(nmf.results,nmf.fit,w,h,feature_score,file="dir/results/NMF_results_rank_3_nrun300_1711_686_bootstrapped_genes_rlm_r.Rdata")

end_time <- Sys.time()
elapsed_time <- end_time - start_time 
elapsed_time
