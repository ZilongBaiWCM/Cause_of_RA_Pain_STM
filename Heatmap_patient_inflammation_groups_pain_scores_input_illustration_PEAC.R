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
library(pheatmap)

setwd("/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/")

peac_dir <- '/Users/baizilong/Documents/Dana-dataset/Validation/'
# Loading PEAC meta data
meta_fn <- 'PEAC FASTQ METADATA.xlsx'
meta_data <- read_excel(paste(peac_dir, meta_fn, sep=""), sheet = 2)
# Focus on "fibroid"

vas_scores <- meta_data[c('Technology','Characteristics[VAS]', 'Characteristics[pathotype]')]

colnames(vas_scores) <- c('Technology', 'VAS', 'Pathotype')

pain_scores <- vas_scores

pain_scores['Inflammation Level'] <- 'NaN'

pain_scores[pain_scores$Pathotype=='lymphoid', 'Inflammation Level'] <- '1.High'
pain_scores[pain_scores$Pathotype=='myeloid', 'Inflammation Level'] <- '2.Mixed'
pain_scores[pain_scores$Pathotype=='fibroid' | pain_scores$Pathotype=='ungraded', 'Inflammation Level'] <- '3.Low'

pain_scores <- unique(pain_scores)

pain_scores_infla <- pain_scores

pain_scores_infla_sorted <- pain_scores_infla %>% arrange(`Inflammation Level`, VAS)

# Load gene expression (bulk, all) from PEAC
exp_fn <- 'batch_corrected_unfiltered_with_gene_log_cpm_peac_04192022.txt'
exp_data <- read.table(paste(peac_dir, exp_fn, sep=""), header=T, row.names =1)

row.names(exp_data) <- gsub("\\..*","",exp_data$gene)

# Load 2227 gene list from HSS.
data_dir = "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/"
svfn = paste(data_dir, 'Generated_Gene_Lists/bulk_padj_<0.01log2FC_<0.csv', sep="")
gene_list <- read.csv(svfn)
gene_list <- gene_list$X

gex_sub <- exp_data[row.names(exp_data) %in% gene_list, pain_scores_infla_sorted$Technology]

# sort gex_sub genes/rows with ascending Laplacian Scores.
# ---- Load Laplacian Scores ...
ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores', filter_condition,'.csv', sep=""))
row.names(ls_result) <- ls_result$X

ls_result <- ls_result[row.names(ls_result) %in% row.names(gex_sub),]

attach(ls_result)
ls_sorted <- ls_result[order(Laplacian_Score),]
detach(ls_result)

gex_sorted <- gex_sub[row.names(ls_sorted), pain_scores_infla_sorted$Technology]

#colnames(pain_scores_infla_sorted) <- c('Pain Score', 'Inflammation Level')
pain_scores_infla_sorted <- data.frame(pain_scores_infla_sorted)
row.names(pain_scores_infla_sorted) <- pain_scores_infla_sorted$Technology
pain_scores_infla_sorted <- pain_scores_infla_sorted %>% rename('Inflammation Level' = 'Inflammation.Level')
annotation_col <- pain_scores_infla_sorted[c('VAS', 'Inflammation Level')]
annotation_col <- data.frame(annotation_col)
annotation_col <- annotation_col %>% rename('Inflammation Level' = 'Inflammation.Level')

top_s <- 738
ls_sorted["Pain Associated"] <- "Less"
ls_sorted[1:top_s, "Pain Associated"] <- "Yes"

gex_sorted <- gex_sorted[row.names(ls_sorted),]

annotation_row <- data.frame(`Laplacian Score` = ls_sorted$Laplacian_Score,
                             `Pain Associated` = ls_sorted$`Pain Associated`)
row.names(annotation_row) <- row.names(ls_sorted)
colnames(annotation_row) <- c('Laplacian Score', 'Pain Associated')

annotation_colors = list(
  `Inflammation Level` = c(`1.High`="orange", `2.Mixed`="brown", `3.Low`="green"),
  `VAS` = colorRampPalette(c("yellow","white","purple"))(60),
  `Laplacian Score` = colorRampPalette(c("gold", "yellow", "white","blue"))(100),
  `Pain Associated` = c(Yes="pink", Less="green"))

save_plot_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper/'
fn <- "Heatmap.low_inflammatory_genes_patients_pain_scores_inflammation_levels_in_PEAC_Fig2G.png"
pheatmap(gex_sorted, color=colorRampPalette(c("blue", "white","red"))(60), scale="row", 
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         fontsize_col = 5,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         treeheight_row = FALSE,
         #main="G",
         #file=paste(save_plot_dir,fn,sep="")
         )

#gex_ENSG_Symbols <- gex[,c('ENSG','Symbol')]
#write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols.csv')
