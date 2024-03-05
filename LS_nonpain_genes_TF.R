# This script aims to find pathways that are "similar" to the Neuroactive Ligand-Receptor Interaction.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
rm(list=ls())
library(magrittr)
library(dplyr)
library(kim)
library(pheatmap)
library(RColorBrewer)
#BiocManager::install("M3C")

#library(M3C)

#help(clustersim)
#help(M3C)

input_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Generated_Gene_lists/"

df_genes <- read.csv(paste(input_dir, 'bulk_padj_<0.01log2FC_<0.csv',sep=""))
genes_lowinf <- df_genes$X
#rm(df_genes)

# Load genes collectively significantly correlate with pain scores in bulk RNA-seq data.
ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
s <- 815
file<-'Laplacian_Scores_of_zscores_bulkonly_top815_in_sorted_769.csv'
laplacian_score_top <- read.csv(paste(ls_dir,file,sep=""))
row.names(laplacian_score_top) <- laplacian_score_top$X

laplacian_score_top <- laplacian_score_top[,2:4]

data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"
gex_ENSG_Symbols <- read.csv(paste(data_dir, 'gex_ENSG_Symbols.csv', sep=""))
row.names(gex_ENSG_Symbols) <- gex_ENSG_Symbols$X

# Convert row names to a new column
#sc_df <- sc_df %>% rownames_to_column()
#gex_ENSG_Symbols <- gex_ENSG_Symbols %>% rownames_to_column()
