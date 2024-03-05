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

# Load the required library
library(base)
#library(h5)
library(rhdf5)


# Specify the path to your .rds.gz file
file_path <- paste(data_dir, 'GSE201546_gpallcelltypes', sep='')

# Read the RDS file
hDRG_scdata <- readRDS(file_path)

file_path <- paste(data_dir, 'GSM6067401_SAM24392087_filtered_feature_bc_matrix.h5', sep='')

hDRG_sample = h5read(file_path, 'matrix') #h5ls(file_path)

# Load human specimens
library(readxl)
metatable <- read_xlsx(paste(data_dir, "41467_2023_36014_MOESM3_ESM.xlsx", sep=''))
unique(metatable$Species)

humans <- metatable[metatable$Species == "Homo sapiens", ]
human_samples = unique(c(humans$SampleID))

print(paste('Number of nuclei from human specimens:', sum(hDRG_scdata$orig.ident %in% human_samples)))

h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)

samples_list <- c()
for (i in 1:length(h5_files)){
  print(h5_files[i])
  tmp <- strsplit(h5_files[i], "_")
  samples_list[length(samples_list) + 1]= tmp[[1]][5]
}
samples_list = unique(samples_list)
print(samples_list)

intersect(unique(hDRG_scdata$orig.ident), unique(metatable$SampleID))

intersect(samples_list, human_samples)
