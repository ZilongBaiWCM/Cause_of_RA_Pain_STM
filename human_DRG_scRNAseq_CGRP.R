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
#colnames(umi) <- paste(human_meta$celltype, human_nuclei,sep="_")

calca_expression <- umi['Calca', ]

human_meta$calca_expression <- calca_expression[rownames(human_meta)]

human_meta$celltype <- as.factor(human_meta$celltype)

num_groups = length(unique(human_meta$celltype))
rainbow_colors <- rainbow(num_groups)

# Create a boxplot with automatically set colors
CALCA_in_sorted <-  ggplot(human_meta, aes(x = celltype, y = calca_expression, fill = celltype)) +
  geom_boxplot(outlier.shape = NA) + #(outlier.shape = NA) +
  geom_jitter(aes(color = celltype), width = 0.2, size=1) +
  #coord_cartesian(ylim = c(0, 6)) +
  scale_fill_manual(values = rainbow_colors) +
  theme(#axis.title.x = element_blank(),  # Disable x-axis title
        #axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text = element_text(size = 7),
        legend.position = "none") +
  labs(x = "Cell Type", y = "Transcripts Per Million (TPMs)")

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

#ggsave(filename = paste(save_dir, 'CALCA_hDRG.pdf', sep=""), plot = CALCA_in_sorted, dpi = 300 ) #, width = 2.2, height = 2.9, dpi = 180, units = "in")

library(svglite)
ggsave(filename = paste(save_dir, 'CALCA_hDRG.svg', sep=""), plot = CALCA_in_sorted, device="svg", width = 7.3, height = 5) #, width = 2.2, height = 2.9, dpi = 180, units = "in")


