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

data <- read_excel('Pain Scores for ZIlong 5.18.23.xls')

data <- data.frame(data)

row.names(data) <- data$Study.ID

VAS_pain <- data$How.much.pain.have.you.felt.due.to.your.rheumatoid.arthritis.during.the.last.week.
  
patient_IDs_low_infla <- read.csv('low_inflammatory_patient_study_ids.csv')
patient_IDs_high_infla <- read.csv('high_inflammatory_patient_study_ids.csv')
colnames(patient_IDs_high_infla) <- colnames(patient_IDs_low_infla)

patient_IDs <- rbind(patient_IDs_low_infla, patient_IDs_high_infla)

meta_pain_scores <- data[data$Study.ID %in% patient_IDs$Study.IDs,c('Study.ID','How.much.pain.have.you.felt.due.to.your.rheumatoid.arthritis.during.the.last.week.')]

pain_scores <- meta_pain_scores$How.much.pain.have.you.felt.due.to.your.rheumatoid.arthritis.during.the.last.week.

names(pain_scores) <- paste('RA',row.names(meta_pain_scores),'.SYN',sep="")
pain_scores <- data.frame(pain_scores)

patient_IDs_low_infla_v <- paste('RA',patient_IDs_low_infla$Study.IDs, '.SYN', sep="")
patient_IDs_high_infla_v <- paste('RA',patient_IDs_high_infla$Study.IDs, '.SYN', sep="")

pain_scores['inflammation_level'] <- NaN

pain_scores[row.names(pain_scores) %in% patient_IDs_low_infla_v, 'inflammation_level'] <- 'Low (combined)'
pain_scores[row.names(pain_scores) %in% patient_IDs_high_infla_v, 'inflammation_level'] <- 'High'

pain_scores_infla <- pain_scores

pain_scores_infla_sorted <- pain_scores_infla %>% arrange(inflammation_level, pain_scores)

# Load gene expression (bulk, all)
gex <- read.delim('dat.DESeq_norm.combat.45samples.with_gene_names.txt', header = TRUE, sep = "\t")

gex <- separate(gex, Gene.ID, into = c("ENSG", "Symbol"), sep = " (?=[^ ]+$)")

row.names(gex) <- gsub("\\..*","",gex$ENSG)

# Load gene list filtered by DEA
data_dir = "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/"
svfn = paste(data_dir, 'Generated_Gene_Lists/bulk_padj_<0.01log2FC_<0.csv', sep="")
gene_list <- read.csv(svfn)
gene_list <- gene_list$X

gex_sub <- gex[gene_list, row.names(pain_scores_infla_sorted)]

# sort gex_sub genes/rows with ascending Laplacian Scores.
# ---- Load Laplacian Scores ...
ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores', filter_condition,'.csv', sep=""))
row.names(ls_result) <- ls_result$X

attach(ls_result)
ls_sorted <- ls_result[order(Laplacian_Score),]
detach(ls_result)

gex_sorted <- gex_sub[row.names(ls_sorted), row.names(pain_scores_infla_sorted)]

colnames(pain_scores_infla_sorted) <- c('Pain Score', 'Inflammation Level')
annotation_col <- pain_scores_infla_sorted

top_s <- 815
ls_sorted["Selected"] <- "No"
ls_sorted[1:top_s, "Selected"] <- "Yes"

annotation_row <- data.frame(Laplacian_score = ls_sorted$Laplacian_Score,
                             We_Select = ls_sorted$Selected)
row.names(annotation_row) <- row.names(ls_sorted)
colnames(annotation_row) <- c('Laplacian Score', 'GbGMI-identified')

annotation_colors = list(
  `Inflammation Level` = c(High="orange", `Low (combined)`="green"),
  `Pain Score` = colorRampPalette(c("yellow","white","purple"))(60),
  #Laplacian_score = rainbow(100))
  #`Laplacian Score` = colorRampPalette(c("gold", "yellow", "white","blue"))(100),
  `Laplacian Score` = colorRampPalette(c("#B22222", "brown", "white","cyan"))(100),
  `GbGMI-identified` = c(Yes="#B22222", No="cyan"))

save_plot_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper/'
fn <- "Heatmap.low_inflammatory_genes_patients_pain_scores_inflammation_levels_ls_sorted_input_Fig2A.png"
pheatmap(gex_sorted, color=colorRampPalette(c("blue", "white","red"))(60), scale="row", 
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE, 
         show_rownames = FALSE,
         cluster_rows = FALSE,
         treeheight_row = FALSE,
         cellheight=.20, cellwidth = 10,
         #main="G",
         #file=paste(save_plot_dir,fn,sep="")
         )

pain_genes = row.names(annotation_row[annotation_row['GbGMI-identified']=='Yes',])

high_inf = row.names(annotation_col[annotation_col['Inflammation Level']=='High',])

low_inf = row.names(annotation_col[annotation_col['Inflammation Level']=='Low (combined)',])

# Extract GEX for all HSS low and high inflammation with VAS pain scores.
## Z-score
gex_highlow_pain <- gex_sorted[pain_genes, colnames(gex_sorted)]
gex_highlow_pain_zscore <- apply(gex_highlow_pain, 1, scale)
gex_highlow_pain_zscore <- t(gex_highlow_pain_zscore)
colnames(gex_highlow_pain_zscore) <- colnames(gex_highlow_pain)
rownames(gex_highlow_pain_zscore) <- rownames(gex_highlow_pain)
write.csv(gex_highlow_pain_zscore,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/GEX_zscores_bulk_pain815_all_HSSVAS.csv',sep=''))

write.csv(pain_scores,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_vashss_all.csv',sep=''))

# Extract GEX for HSS high inflammation with VAS pain scores.
## Z-score
gex_highinf_pain <- gex_sorted[pain_genes, high_inf]
gex_highinf_pain_zscore <- apply(gex_highinf_pain, 1, scale)
gex_highinf_pain_zscore <- t(gex_highinf_pain_zscore)
colnames(gex_highinf_pain_zscore) <- colnames(gex_highinf_pain)
rownames(gex_highinf_pain_zscore) <- rownames(gex_highinf_pain)
write.csv(gex_highinf_pain_zscore,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/GEX_zscores_bulk_pain815_highinf_HSSVAS.csv',sep=''))

write.csv(pain_scores[high_inf,],paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_vashss_highinf.csv',sep=''))

# Extract GEX for HSS low inflammation with VAS pain scores.
## Z-score
gex_lowinf_pain <- gex_sorted[pain_genes, low_inf]
gex_lowinf_pain_zscore <- apply(gex_lowinf_pain, 1, scale)
gex_lowinf_pain_zscore <- t(gex_lowinf_pain_zscore)
colnames(gex_lowinf_pain_zscore) <- colnames(gex_lowinf_pain)
rownames(gex_lowinf_pain_zscore) <- rownames(gex_lowinf_pain)
write.csv(gex_lowinf_pain_zscore,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/GEX_zscores_bulk_pain815_lowinf_HSSVAS.csv',sep=''))

write.csv(pain_scores[low_inf,],paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_vashss_lowinf.csv',sep=''))

