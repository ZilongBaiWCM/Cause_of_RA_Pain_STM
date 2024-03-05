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

pain_scores[row.names(pain_scores) %in% patient_IDs_low_infla_v, 'inflammation_level'] <- 'Low (combined)'
pain_scores[row.names(pain_scores) %in% patient_IDs_high_infla_v, 'inflammation_level'] <- 'High'

pain_scores_infla <- pain_scores

pain_scores_infla_sorted <- pain_scores_infla %>% arrange(inflammation_level, pain_scores)

# Load gene expression (bulk, all)
gex <- read.delim('dat.DESeq_norm.combat.45samples.with_gene_names.txt', header = TRUE, sep = "\t")

gex <- separate(gex, Gene.ID, into = c("ENSG", "Symbol"), sep = " (?=[^ ]+$)")

row.names(gex) <- gsub("\\..*","",gex$ENSG)

selected_gene = 'TF'

high_inf = row.names(pain_scores[pain_scores['inflammation_level']=='High',])

low_inf = row.names(pain_scores[pain_scores['inflammation_level']=='Low (combined)',])

# Extract GEX for all HSS low and high inflammation with HOOS/KOOS pain scores.
## Z-score
gex_highlow_pain <- gex[gex$Symbol == selected_gene, row.names(pain_scores)]
rownames(gex_highlow_pain) <- c('TF')
write.csv(gex_highlow_pain,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/GEX_TF_bulk.csv',sep=''))

#write.csv(pain_scores[colnames(gex_highlow_pain_zscore),],paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_hooskoos.csv',sep=''))

# Extract GEX for HSS high inflammation with HOOS/KOOS pain scores.
## Z-score
gex_high_pain <- gex[gex$Symbol == selected_gene, high_inf]
rownames(gex_high_pain) <- c('TF')
write.csv(gex_high_pain,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/GEX_TF_bulk_highinf.csv',sep=''))

#write.csv(pain_scores[high_inf,],paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_hooskoos_highinf.csv',sep=''))

# Extract GEX for HSS low inflammation with HOOS/KOOS pain scores.
## Z-score
gex_low_pain <- gex[gex$Symbol == selected_gene, low_inf]
rownames(gex_low_pain) <- c('TF')
write.csv(gex_low_pain,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/GEX_TF_bulk_lowinf.csv',sep=''))

#write.csv(pain_scores[low_inf,],paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_hooskoos_lowinf.csv',sep=''))

