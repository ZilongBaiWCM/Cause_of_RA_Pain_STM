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


excel_file_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'
meta <- read_excel(paste(excel_file_dir, 'Supplemental Table Metadata dat.DESeq_norm.combat.45samples..forpain..xlsx', sep=""))

meta_selected <- meta[meta$`Study ID` %in% patient_IDs$Study.IDs,]

file_name <- 'Data file S4 Patient demographic data and RA disease activity scores and treatments.xlsx'
file_path <- paste(excel_file_dir,file_name, sep="")

my_data <- data.frame(meta_selected, row.names = NULL)
write.xlsx(my_data, file_path, sheetName = "Metadata of RA patients")
