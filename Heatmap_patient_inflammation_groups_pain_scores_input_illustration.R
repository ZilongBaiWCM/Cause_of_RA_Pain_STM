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

data <- read_csv('rheumatoid_arthritis_redcap_raw_data.csv')

data <- data.frame(data)

row.names(data) <- data$study_id

right_hip <- data[data$op_joint == 1,]
left_hip <- data[data$op_joint == 2,]
both_hips <- data[data$op_joint == 3,]
right_knee <- data[data$op_joint == 4,]
left_knee <- data[data$op_joint == 5,]
both_knees <- data[data$op_joint == 6,]
hips_and_knees <- rbind(right_hip, left_hip, both_hips, right_knee, left_knee, both_knees)

koos_pain <- hips_and_knees[hips_and_knees$sc_koos_pn != 999,]
hoos_pain <- hips_and_knees[hips_and_knees$sc_hoos_pn != 999,]

patient_IDs_low_infla <- read.csv('low_inflammatory_patient_study_ids.csv')
patient_IDs_high_infla <- read.csv('high_inflammatory_patient_study_ids.csv')
colnames(patient_IDs_high_infla) <- colnames(patient_IDs_low_infla)

patient_IDs <- rbind(patient_IDs_low_infla, patient_IDs_high_infla)

meta_pain_scores <- data[data$study_id %in% patient_IDs$Study.IDs,c('sc_hoos_pn','sc_koos_pn')]

pain_scores <- c()

for( i in 1:dim(meta_pain_scores)[1]){
  pain_scores <- append(pain_scores, min(meta_pain_scores[i,]))
}

names(pain_scores) <- paste('RA',row.names(meta_pain_scores),'.SYN',sep="")
pain_scores <- data.frame(pain_scores)

patient_IDs_low_infla_v <- paste('RA',patient_IDs_low_infla$Study.IDs, '.SYN', sep="")
patient_IDs_high_infla_v <- paste('RA',patient_IDs_high_infla$Study.IDs, '.SYN', sep="")

pain_scores['inflammation_level'] <- NaN

pain_scores[row.names(pain_scores) %in% patient_IDs_low_infla_v, 'inflammation_level'] <- 'Relatively Low'
pain_scores[row.names(pain_scores) %in% patient_IDs_high_infla_v, 'inflammation_level'] <- 'High'

pain_scores_infla <- pain_scores

pain_scores_infla_sorted <- pain_scores_infla %>% arrange(inflammation_level, pain_scores)

#write.csv(pain_scores,'pain_scores.csv')

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

annotation_colors = list(
  `Inflammation Level` = c(`Relatively Low`="pink", High="#B22222"),
  `Pain Score` = colorRampPalette(c("yellow","white","purple"))(60),
  #Laplacian_score = rainbow(100))
  Laplacian_score = colorRampPalette(c("gold", "yellow", "white","blue"))(100),
  We_Select = c(Yes="gold", No="blue"))

save_plot_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper/'
fn <- "Heatmap.low_inflammatory_genes_patients_pain_scores_inflammation_levels_input_SupFig1B.png"
pheatmap(gex_sorted, color=colorRampPalette(c("green", "white","red"))(60), scale="row", 
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE, show_rownames = FALSE,
         treeheight_row = FALSE,
         #main="G",
         file=paste(save_plot_dir,fn,sep=""))

gex_ENSG_Symbols <- gex[,c('ENSG','Symbol')]
#write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols.csv')
pheatmap(gex_sorted, color=colorRampPalette(c("green", "white","red"))(60), scale="row", 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE, show_rownames = FALSE,
         treeheight_row = FALSE)
