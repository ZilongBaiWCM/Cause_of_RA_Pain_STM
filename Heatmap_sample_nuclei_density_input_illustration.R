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

top_k = 5000

nuc_density <- read.csv(paste('nuc_density_scores_all_patients_',top_k,'_genes.csv'), row.names=1)

meta_dir <- '/Users/baizilong/Documents/Dana-dataset/Bulk Synovial Tissue RNAseq/ResultFiles/RNA_sequencing_result/'

meta_rna <- read.xlsx(paste(meta_dir, '45samples.Metadata..xls', sep=""), sheetIndex = 1)

row.names(meta_rna) <- meta_rna$Study.ID

print(dim(meta_rna))
meta_rna <- meta_rna[meta_rna$actual.clusters == 1 | meta_rna$actual.clusters == 2 | meta_rna$actual.clusters == 3, ] # All clusters taken into consideration. Should not have changed the dimensionality of meta_rna.
print(dim(meta_rna))

meta_rna$Study.ID <- paste('RA',meta_rna$Study.ID,'.SYN',sep="")

row.names(meta_rna) <- meta_rna$Study.ID

study.ids <- intersect(row.names(meta_rna), row.names(nuc_density))

meta_rna <- meta_rna[study.ids, ]
nuc_density <- nuc_density[study.ids, ]
names(nuc_density) <- study.ids

nuc_density <- data.frame(nuc_density)
colnames(nuc_density) <- c("nuc_density_scores")

nuc_density['cluster.id'] <- meta_rna$actual.clusters

nuc_density['Inflammation Level'] <- 'NaN'

nuc_density[nuc_density$cluster.id==1, 'Inflammation Level'] <- 'Low'
nuc_density[nuc_density$cluster.id==2, 'Inflammation Level'] <- 'Mixed'
nuc_density[nuc_density$cluster.id==3, 'Inflammation Level'] <- 'High'

nuc_density <- unique(nuc_density)

nuc_density_infla <- nuc_density

#nuc_density_infla_sorted <- nuc_density_infla %>% arrange(cluster.id, nuc_density_scores)
nuc_density_infla_sorted <- nuc_density_infla %>% arrange(nuc_density_scores)

# Load gene expression (bulk, top 5000) from HSS
exp_data <- read.csv(paste('GEX_bulk_vs_nuc_all_patients_top_',top_k,'_genes.csv',sep=""), row.names =1)

gex_sub <- exp_data[row.names(exp_data), row.names(nuc_density_infla_sorted)]

gex_sorted <- gex_sub

# nuc_density_infla_sorted<- data.frame(nuc_density_infla_sorted)
# nuc_density_infla_sorted <- nuc_density_infla_sorted %>% rename('Inflammation Level' = 'Inflammation.Level')

annotation_col <- nuc_density_infla_sorted[c('nuc_density_scores')]
annotation_col <- data.frame(annotation_col)
colnames(annotation_col) <- c('Nuclei Density')

annotation_colors = list(
  #`Inflammation Level` = c('Low'="green", 'Medium'="yellow", 'High'="purple"),
  `Nuclei Density` = colorRampPalette(c("yellow","white","purple"))(60))

save_plot_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper/'
fn <- "Heatmap.low_inflammatory_genes_patients_nuclei_density_scores_in_HSS_SupFig2A.png"
pheatmap(gex_sorted, color=colorRampPalette(c("blue", "white","red"))(60), scale="row", 
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         fontsize_col = 8,
         #cluster_rows = FALSE,
         show_rownames = FALSE,
         treeheight_row = FALSE,
         cellheight=.10, cellwidth = 15,
         #main="G",
         file=paste(save_plot_dir,fn,sep="")
         )

#gex_ENSG_Symbols <- gex[,c('ENSG','Symbol')]
#write.csv(gex_ENSG_Symbols, 'gex_ENSG_Symbols.csv')
