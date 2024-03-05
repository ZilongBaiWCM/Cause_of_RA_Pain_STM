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

setwd("/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/")

peac_dir <- '/Users/baizilong/Documents/Dana-dataset/Validation/'
# Loading PEAC meta data
meta_fn <- 'PEAC FASTQ METADATA.xlsx'
meta_data <- read_excel(paste(peac_dir, meta_fn, sep=""), sheet = 2)
# Focus on "fibroid"

vas_scores <- meta_data[c('Technology','Characteristics[VAS]', 'Characteristics[pathotype]')]

colnames(vas_scores) <- c('Technology', 'VAS', 'Pathotype')

pain_scores <- vas_scores

pain_scores <- data.frame(pain_scores[!duplicated(pain_scores), ])

rownames(pain_scores) <- pain_scores$Technology

# Load gene expression (bulk, all)
exp_fn <- 'batch_corrected_unfiltered_with_gene_log_cpm_peac_04192022.txt'
gex <- read.table(paste(peac_dir, exp_fn, sep=""), header=T, row.names =1)

row.names(gex) <- gsub("\\..*","",gex$gene)

# How do we get the input to the algorithm? How to have a fair comparison? Should we get the DEA specifically from this dataset?
# Use the previously identified low-inflammatory genes.
# Load 2227 gene list from HSS.
data_dir = "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/"
svfn = paste(data_dir, 'Generated_Gene_Lists/bulk_padj_<0.01log2FC_<0.csv', sep="")
gene_list <- read.csv(svfn)
gene_list <- gene_list$X

gex_sub <- gex[row.names(gex) %in% gene_list, row.names(pain_scores)]

gex_zscores <- apply(gex_sub, MARGIN = 1, FUN = scale )
gex_zscores <- t(gex_zscores)

colnames(gex_zscores) <- colnames(gex_sub)

# Pre-processing to get PEAC patients with VAS pain scores in lymphoid (high inflammation).
pain_scores_lymphoid <- pain_scores[pain_scores$Pathotype=='lymphoid', ]

#write.csv(pain_scores_lymphoid,'pain_scores_peac_lymphoid.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

# Gene expression must be gene-wise z-scored over all the investigated patients within a given group.

gex_zscores_lymphoid <- gex_zscores[, rownames(pain_scores_lymphoid)]

write.csv(gex_zscores_lymphoid, 'GEX_zscores_bulk_lymphoid_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')

###################
###################

# Pre-processing to get PEAC patients with VAS pain scores in myeloid (mixed inflammation).
pain_scores_myeloid <- pain_scores[pain_scores$Pathotype=='myeloid', ]

gex_zscores_myeloid <- gex_zscores[, rownames(pain_scores_myeloid)]

gex_zscores_myeloid <- data.frame(gex_zscores_myeloid)

write.csv(gex_zscores_myeloid, 'GEX_zscores_bulk_myeloid_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')

###################
###################

# Pre-processing to get PEAC patients with VAS pain scores in fibroid/ungraded (low inflammation).
pain_scores_fibroidungraded <- pain_scores[pain_scores$Pathotype=='fibroid' | pain_scores$Pathotype=='ungraded', ]

#write.csv(pain_scores_fibroidungraded,'pain_scores_peac_fibroidungraded.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

gex_zscores_fibroidungraded <- gex_zscores[, rownames(pain_scores_fibroidungraded)]

gex_zscores_fibroidungraded <- data.frame(gex_zscores_fibroidungraded)

write.csv(gex_zscores_fibroidungraded, 'GEX_zscores_bulk_fibroidungraded_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')

###################
###################

# Pre-processing to get PEAC patients with VAS pain scores in ALL (lymphoid, myeloiod, fibroid, and ungraded)..
pain_scores_all <- pain_scores[, ] # All patients

#write.csv(pain_scores_all,'pain_scores_peac_all.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

gex_zscores_all <- gex_zscores[, rownames(pain_scores_all)]

gex_zscores_all <- data.frame(gex_zscores_all)

write.csv(gex_zscores_all, 'GEX_zscores_bulk_all_PEAC_v2.csv')

###### 
# TODO: Myeloid + Fibroid
# Pre-processing to get PEAC patients with VAS pain scores in myeloid/fibroid (low inflammation).
pain_scores_myefib <- pain_scores[pain_scores$Pathotype=='myeloid' | pain_scores$Pathotype=='fibroid', ]

#write.csv(pain_scores_myefib,'pain_scores_peac_myefib.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

gex_zscores_myefib <- gex_zscores[, rownames(pain_scores_myefib)]

gex_zscores_myefib <- data.frame(gex_zscores_myefib)

write.csv(gex_zscores_myefib, 'GEX_zscores_bulk_myefib_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')

######
# TODO: Fibroid alone, without ungraded
# Pre-processing to get PEAC patients with VAS pain scores in fibroid (low inflammation).
pain_scores_fibroid <- pain_scores[pain_scores$Pathotype=='fibroid', ]

#write.csv(pain_scores_fibroid,'pain_scores_peac_fibroid.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

gex_zscores_fibroid <- gex_zscores[, rownames(pain_scores_fibroid)]

gex_zscores_fibroid <- data.frame(gex_zscores_fibroid)

write.csv(gex_zscores_fibroid, 'GEX_zscores_bulk_fibroid_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')

######
# TODO: Myeloid + Fibroid + Ungraded
# Pre-processing to get PEAC patients with VAS pain scores in myeloid/fibroid/ungraded (low inflammation).
pain_scores_myefibun <- pain_scores[pain_scores$Pathotype=='myeloid' | pain_scores$Pathotype=='fibroid' | pain_scores$Pathotype=='ungraded', ]

#write.csv(pain_scores_myefibun,'pain_scores_peac_myefibun.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

gex_zscores_myefibun <- gex_zscores[, rownames(pain_scores_myefibun)]

gex_zscores_myefibun <- data.frame(gex_zscores_myefibun)

write.csv(gex_zscores_myefibun, 'GEX_zscores_bulk_myefibun_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')

######
# TODO: Lymphoid + Fibroid + Ungraded
# Pre-processing to get PEAC patients with VAS pain scores in lymphoid/fibroid/ungraded (all inflammation presented).
pain_scores_lymfibun <- pain_scores[pain_scores$Pathotype=='lymphoid' | pain_scores$Pathotype=='fibroid' | pain_scores$Pathotype=='ungraded', ]

write.csv(pain_scores_lymfibun,'pain_scores_peac_lymfibun.csv') # This file will be used to generate the pair-wise pain-based similarity graph.

gex_zscores_lymfibun <- gex_zscores[, rownames(pain_scores_lymfibun)]

gex_zscores_lymfibun <- data.frame(gex_zscores_lymfibun)

write.csv(gex_zscores_lymfibun, 'GEX_zscores_bulk_lymfibun_PEAC_v2.csv')

# gex_ENSG_Symbols <- gex[,c('gene','SYMBOL.1')]
# write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols_PEAC.csv')
