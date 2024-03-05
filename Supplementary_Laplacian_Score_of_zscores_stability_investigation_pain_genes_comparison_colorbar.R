# This script aims to find pathways that are "similar" to the Neuroactive Ligand-Receptor Interaction.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
rm(list=ls())

#library(SingleCellExperiment)
library(ggplot2)
library(magrittr)
library(DESeq2)
library(dplyr)
#install.packages("factoextra")
library(factoextra)
#BiocManager::install("clusterProfiler")
library(enrichplot)
library(org.Hs.eg.db) # Gene Ontology
library(clusterProfiler)
#BiocManager::install('M3C')
library(M3C) # TSNE
library(umap) # umap
library(cluster)
#install.packages('amap')
library(amap) # Kmeans
library(tidyverse)  # data manipulation
library(dendextend) # for comparing two dendrograms
# load library
library(ggrepel) # Volcano plot
library(tidyverse)
library(readxl)
#install.packages('xlsx')
library(xlsx)

library(ggplot2)

library(Rtsne)

library(pheatmap)

data_dir = '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/'
df_ps <- read.csv(paste(data_dir, 'pain_scores.csv', sep=""), row.names=1)
N = dim(df_ps)[1]

gex_zscores <- read.csv(paste(data_dir, 'GEX_zscores_bulk_padj_<0.01log2FC_<0.csv',sep=""), row.names=1) # TODO: we need to get different z-scored gene expression as the patient group is changed.

genes_background <- row.names(gex_zscores)

input_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'

settings <- c(paste('Exclude', row.names(df_ps)), "22 rel low-inf") # The last row/column is for the original patient list.

mat_jaccard <- matrix(, nrow = N+1, ncol = N+1)
row.names(mat_jaccard) <- settings
colnames(mat_jaccard) <- settings
mat_fishers <- matrix(, nrow = N+1, ncol = N+1)
row.names(mat_fishers) <- settings
colnames(mat_fishers) <- settings

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# ONGOING. Incomplete.
perp <- 6 # Here the perplexity used for LS computation with one patient excluded.
for (i in 1:N){
  df_pg_noi <- read.csv(paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/pain_related_genes_GbGMI_perplexity',perp,'_exclude_patient_',i,'.csv',sep=''), row.names=1)
  genes_noi <- row.names(df_pg_noi)
  for (j in 1:N){
    df_pg_noj <- read.csv(paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/pain_related_genes_GbGMI_perplexity',perp,'_exclude_patient_',j,'.csv',sep='') , row.names=1)
    genes_noj <- row.names(df_pg_noj)
    
    mat_jaccard[i, j] <- jaccard(genes_noi, genes_noj)
    
    dat <- data.frame(
      "noj_yes" = c(length(intersect(genes_noj, genes_noi)), 
                           length(intersect(genes_noj, setdiff(genes_background, genes_noi)))),
      "noj_no" = c(length(intersect(setdiff(genes_background, genes_noj), genes_noi)), 
                          length(intersect(setdiff(genes_background, genes_noj), setdiff(genes_background, genes_noi)))),
      row.names = c("noi_yes", "noi_no"),
      stringsAsFactors = FALSE
    )
    colnames(dat) <- c("noj_yes", "noj_no")
    test <- fisher.test(dat)
    test$p.value
    
    mat_fishers[i, j] <- test$p.value
  }
}

ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores_of_zscores_bulkonly_top815.csv', sep=""), row.names=1)
genes_all_patients <- row.names(ls_result)

for (i in 1:N){
  df_pg_noi <- read.csv(paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/pain_related_genes_GbGMI_perplexity',perp,'_exclude_patient_',i,'.csv',sep=''), row.names=1)
  genes_noi <- row.names(df_pg_noi)
  
  mat_jaccard[i, N+1] <- jaccard(genes_noi, genes_all_patients)
  mat_jaccard[N+1, i] <- jaccard(genes_all_patients, genes_noi)
  
  dat <- data.frame(
    "noi_yes" = c(length(intersect(genes_noi, genes_all_patients)), 
                  length(intersect(genes_noi, setdiff(genes_background, genes_all_patients)))),
    "noi_no" = c(length(intersect(setdiff(genes_background, genes_noi), genes_all_patients)), 
                 length(intersect(setdiff(genes_background, genes_noi), setdiff(genes_background, genes_all_patients)))),
    row.names = c("ap_yes", "ap_no"),
    stringsAsFactors = FALSE
  )
  colnames(dat) <- c("noi_yes", "noi_no")
  test <- fisher.test(dat)
  test$p.value
  
  mat_fishers[i, N+1] <- test$p.value
  mat_fishers[N+1, i] <- test$p.value

}


mat_jaccard[N+1, N+1] <- jaccard(genes_all_patients, genes_all_patients)

dat <- data.frame(
  "ap_yes" = c(length(intersect(genes_all_patients, genes_all_patients)), 
                length(intersect(genes_all_patients, setdiff(genes_background, genes_all_patients)))),
  "ap_no" = c(length(intersect(setdiff(genes_background, genes_all_patients), genes_all_patients)), 
               length(intersect(setdiff(genes_background, genes_all_patients), setdiff(genes_background, genes_all_patients)))),
  row.names = c("ap_yes", "ap_no"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("ap_yes", "ap_no")
test <- fisher.test(dat)
test$p.value

mat_fishers[N+1, N+1] <- test$p.value

# Present the correlations in heatmap
save_plot_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- "Heatmap.jaccard_similarity_pain_genes_single_patient_excluded_new_colorbar.pdf"
pheatmap(mat_jaccard, #color=colorRampPalette(c("green", "white","red"))(60), scale="row",
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = F,
         cluster_rows = FALSE,
         #treeheight_row = FALSE,
         cellheight=20, cellwidth = 20,
         #main="G",
         file=paste(save_plot_dir,fn,sep=""))

fn <- "Heatmap.fishers_exact_test_pain_genes_single_patient_excluded_new_colorbar.pdf"
pheatmap(mat_fishers, #color=colorRampPalette(c("green", "white","red"))(60), scale="row",
         cluster_cols = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         cluster_rows = FALSE,
         #treeheight_row = FALSE,
         cellheight=20, cellwidth = 20,
         #main="G",
         file=paste(save_plot_dir,fn,sep=""))
# NOTE: All p-values are precisely zero. No heatmap generated for the p-values.

# Present the correlations in boxplot. What is the best way to present correlation coefficients and p-values at the same time???
