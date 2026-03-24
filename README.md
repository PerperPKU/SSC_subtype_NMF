# SSC_subtype_NMF
Using NMF to subtype ASD
## 1. required packages
### data.table (1.17.4)
### MASS (7.3.58.1)
### sfsmisc (1.1.20)
### optparse (1.7.5)
### NMF (0.28)
## 2. install NMF
### install.package('NMF')
Normally, the installation should take no more than 20 minutes.
## 3. Identify genes related traits using robost linear regression
1.rlm_single_task.R
## 4. Run NMF using trait_related gene expression matrix (chooose optimal rank)
2. NMF_random_n30_1711ASD_5pheno_686_genes_overlap_genes.R
## 5. Run NMF using trait_related gene expression matrix using optimal rank
3.run_NMF_with_rank_3_1711ASD_686genes_rlm


