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

row.names(data) <- paste('RA',data$Study.ID,'.SYN',sep="") 

# Compare Low-Inflammation group
ra_vas <- read.csv(
  paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_hss_ra_vas_lowinf.csv',sep=''),
  row.names=1)

condition_vas <- read.csv(
  paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_hss_condition_vas_lowinf.csv',sep=''),
  row.names=1)

hooskoos <- read.csv(
  paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/pain_scores_hooskoos_lowinf.csv',sep=''),
  row.names=1)

patients<- rownames(hooskoos)

hooskoos <- hooskoos[patients,]
ra_vas <- ra_vas[patients,]
condition_vas <- condition_vas[patients,]

data <- data[patients, c('KOOS.Pain', 'HOOS.Pain')]

data_comp <- data.frame(hooskoos = hooskoos$pain_scores, ra_vas = ra_vas$pain_scores, condition_vas = condition_vas$pain_scores)
rownames(data_comp) <- patients

kendall.result_hooskoosVSra <- cor.test(data_comp$hooskoos, data_comp$ra_vas, method="kendall")
kendall.result_hooskoosVScondition <- cor.test(data_comp$hooskoos, data_comp$condition_vas, method="kendall")



# Compare High-Inflammation group

# Compare All