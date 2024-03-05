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
file_path <- paste(data_dir, 'GSE201654_mousegpcynohumanDRG', sep='')

# Read the RDS file
DRG_scdata <- readRDS(file_path)

file_path <- paste(data_dir, 'GSE201654_RAW/','GSM6067401_SAM24392087_filtered_feature_bc_matrix.h5', sep='')

DRG_sample = h5read(file_path, 'matrix') #h5ls(file_path)

# Load human specimens
library(readxl)
metatable <- read_xlsx(paste(data_dir, "41467_2023_36014_MOESM3_ESM.xlsx", sep=''))
unique(metatable$Species)

humans <- metatable[metatable$Species == "Homo sapiens", ]
human_samples = unique(c(humans$SampleID))

print(paste('Number of nuclei from human specimens:', sum(DRG_scdata$orig.ident %in% human_samples)))

h5_files <- list.files(paste(data_dir, 'GSE201654_RAW/', sep=""), pattern = "\\.h5$", full.names = TRUE)

samples_list <- c()
for (i in 1:length(h5_files)){
  print(h5_files[i])
  tmp <- strsplit(h5_files[i], "_")
  samples_list[length(samples_list) + 1]= tmp[[1]][6]
}
samples_list = unique(samples_list)
print(samples_list)

print('Human samples from meta table MOESM3:')
print(human_samples)

print('Human in both data table and meta table MOESM3:')
print(intersect(unique(DRG_scdata$orig.ident), human_samples))

unique_species = unique(DRG_scdata@meta.data$Species)
#"Mouse"     "Human"     "Cyno"      "Guineapig"

human_meta = DRG_scdata@meta.data[DRG_scdata@meta.data$Species == "Human",]

human_nuclei = row.names(human_meta)

normalized_data = DRG_scdata@assays$RNA@data

human_norm_data = normalized_data[, c(human_nuclei)]

# # Convert mouse gene symbols into human gene symbols
# library(nichenetr)
# mgenes <- row.names(human_norm_data)
# hgenes <- mgenes %>% convert_mouse_to_human_symbols()
# names(hgenes) <- mgenes
# na_indices <- which(is.na(hgenes))
# #human_norm_data <- human_norm_data[-names(na_indices), ]
# human_norm_data <- human_norm_data[!(rownames(human_norm_data) %in% names(na_indices)), ]
# #names(na_indices) <- NULL

human_cell_types <- human_meta$celltype
names(human_cell_types) <- row.names(human_meta)

umi <- human_norm_data
colnames(umi) <- paste(human_meta$celltype, human_nuclei,sep="_")

library(Seurat)
umi_seurat <- CreateSeuratObject(counts = umi, project = "normalized")
umi_seurat@meta.data$groups <- as.factor(human_cell_types)

cell_types_levels <- levels(umi_seurat)
de.markers.list <- list()
ctype.de.markers.list <- list()

library(ggrepel)

ls_dir <- "/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/Ligand_Receptor_hDRG_snRNAseq/"

neuron_type_vs_others_combined <- data.frame()

for (c_type in cell_types_levels){
  #for (c_type in cell_types[cell_types != f_type]){
  de.markers <- FindMarkers(umi_seurat, ident.1 = c_type, ident.2 = NULL)
  
  de_fn <- paste(ls_dir, 'Seurat_',c_type,'_vs_AllOthers_hDRG_snRNAseq.csv', sep="")
  
  write.csv(de.markers, de_fn)
  
  attach(de.markers)
  dea_sc_genes_pval <- rownames(de.markers[which(avg_log2FC>0 & p_val < 0.05),])
  dea_sc_genes_padj <- rownames(de.markers[which(avg_log2FC>0 & p_val_adj < 0.05),])
  detach(de.markers)
  
  fn_de_cell_type_genes <- paste(ls_dir, 'Seurat_',c_type,'_vs_AllOthers_padj_genes_hDRG_snRNAseq.csv', sep="")
  
  write.csv(dea_sc_genes_padj, fn_de_cell_type_genes)
  
  print(paste('Current cell type:', c_type))
  
  print(paste('Number of differentially expressed genes in ', c_type))
  print(length(dea_sc_genes_padj))
  
  print(' "GRIK4" %in% dea_sc_genes_padj" ')
  print("Grik4" %in% dea_sc_genes_padj)
  
  print(' "Tgfbr1" %in% dea_sc_genes_padj" ')
  print("Tgfbr1" %in% dea_sc_genes_padj)
  
  print(' "NTN4" %in% dea_sc_genes_padj" ')
  print("Ntn4" %in% dea_sc_genes_padj)
  
  print(' "UNC5A" %in% dea_sc_genes_padj" ')
  print("Unc5a" %in% dea_sc_genes_padj)
  
  neuron_type_vs_others <- data.frame(dea_sc_genes_padj)
  colnames(neuron_type_vs_others) <- 'mgenes'
  neuron_type_vs_others$cell_type <- c_type
  
  neuron_type_vs_others_combined <- rbind(neuron_type_vs_others_combined, neuron_type_vs_others)
}

neuron_type_vs_others_combined$hgenes <- convert_mouse_to_human_symbols(neuron_type_vs_others_combined$mgenes)

fn_de_cell_type_genes_combined <- paste(ls_dir, 'Seurat_neuron_vs_AllOthers_padj_genes_hDRG_snRNAseq.csv', sep="")

write.csv(neuron_type_vs_others_combined, fn_de_cell_type_genes_combined)