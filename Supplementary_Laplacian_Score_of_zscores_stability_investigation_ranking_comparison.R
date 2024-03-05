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

input_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'

settings <- c(paste('Exclude', row.names(df_ps)), "HSS 22 Low-inflammation Patients") # The last row/column is for the original patient list.

mat_estimate <- matrix(, nrow = N+1, ncol = N+1)
row.names(mat_estimate) <- settings
colnames(mat_estimate) <- settings
mat_p.value <- matrix(, nrow = N+1, ncol = N+1)
row.names(mat_p.value) <- settings
colnames(mat_p.value) <- settings

for (i in 1:N){
  df_ls_noi <- read.csv(paste(input_dir,'Laplacian_Scores',filter_condition,'_exclude_col',i-1,'.csv', sep=""), row.names = 1)
  for (j in 1:N){
    df_ls_noj <- read.csv(paste(input_dir,'Laplacian_Scores',filter_condition,'_exclude_col',j-1,'.csv', sep=""), row.names = 1)
    comp <- merge(df_ls_noi, df_ls_noj,
                  by = 'row.names', all = TRUE)
    res<-cor.test(comp$Laplacian_Score.x,comp$Laplacian_Score.y, method="spearman")
    mat_estimate[i, j] <- res$estimate
    mat_p.value[i, j] <- res$p.value
  }
}

df_ls_all <- read.csv(paste(input_dir,'Laplacian_Scores',filter_condition,'.csv', sep=""), row.names = 1)
for (i in 1:N){
  df_ls_noi <- read.csv(paste(input_dir,'Laplacian_Scores',filter_condition,'_exclude_col',i-1,'.csv', sep=""), row.names = 1)
  comp <- merge(df_ls_noi, df_ls_all,
                by = 'row.names', all = TRUE)
  res<-cor.test(comp$Laplacian_Score.x,comp$Laplacian_Score.y, method="spearman")
  mat_estimate[i, N+1] <- res$estimate
  mat_estimate[N+1, i] <- res$estimate
  mat_p.value[i, N+1] <- res$p.value
  mat_p.value[N+1, i] <- res$p.value
}

comp <- merge(df_ls_all, df_ls_all,
              by = 'row.names', all = TRUE)
res<-cor.test(comp$Laplacian_Score.x,comp$Laplacian_Score.y, method="spearman")
mat_estimate[N+1, N+1] <- res$estimate
mat_p.value[N+1, N+1] <- res$p.value

# Present the correlations in heatmap
save_plot_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- "Heatmap.correlation_coefficient_ls_based_gene_prioritizations_single_patient_excluded_new.png"
pheatmap(mat_estimate, #color=colorRampPalette(c("green", "white","red"))(60), scale="row",
         cluster_cols = FALSE,
         show_rownames = F,
         show_colnames = F,
         legend = F,
         cluster_rows = FALSE,
         #treeheight_row = FALSE,
         cellheight=20, cellwidth = 20,
         #main="G",
         file=paste(save_plot_dir,fn,sep=""))

# fn <- "Heatmap.correlation_p.value_ls_based_gene_prioritizations_single_patient_excluded_new.png"
# pheatmap(mat_p.value, #color=colorRampPalette(c("green", "white","red"))(60), scale="row",
#          cluster_cols = FALSE,
#          show_rownames = T,
#          cluster_rows = FALSE,
#          #treeheight_row = FALSE,
#          cellheight=15, cellwidth = 15,
#          #main="G",
#          file=paste(save_plot_dir,fn,sep=""))
# NOTE: All p-values are precisely zero. No heatmap generated for the p-values.

# Present the correlations in boxplot. What is the best way to present correlation coefficients and p-values at the same time???
