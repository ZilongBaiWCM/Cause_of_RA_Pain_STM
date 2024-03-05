# Zilong Bai, PhD, composed this code.
# For paper: Machine Learning Reveals Synovial Fibroblast Genes Associated with Pain Affect Sensory Nerve Growth in Rheumatoid Arthritis
# Zilong Bai is the first author of this paper.
#
rm(list=ls())
#setwd(getwd())
setwd("/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/Fibroblast_subtyping/amp_phase1_ra-master/R")

library(devtools)
library(igraph)
library(pheatmap)
library(Rtsne)
library(viridis)
library(RColorBrewer)
library(loe)
library(limma)
library(pbapply)
library(vegan)
library(cluster)
library(dplyr)
library(ggplot2)
library(CCA)
library(superheat)
#library(biomaRt)
#library(edgeR) # To use the cpm function for computing log2cpm from integer counts.
require(gdata)

source("meta_colors.R")
source("pure_functions.R")


sc_dir <- "/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/Fibroblast_subtyping/amp_phase1_ra-master/data/"
fn_nd = paste(sc_dir,"celseq_synovium_log2_5265cells_paper.rds", sep="")

fn_meta = paste(sc_dir,"celseq_synovium_meta_5265cells_paper.rds", sep="")

umi_5265cells <- readRDS(fn_nd)
  
meta_5265cells <- readRDS(fn_meta)
rownames(meta_5265cells) <- meta_5265cells$cell_name

meta_RA <- meta_5265cells[meta_5265cells$disease=="RA", ]
umi_RA <- umi_5265cells[, rownames(meta_RA)]

umi <- umi_RA # Only focus on RA cells in this analysis.
meta <- meta_RA

# Filter umi Seurat object with genes of interest identified by our algorithm

dir_goi<- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/"

fn_goi <- "Laplacian_Scores_of_zscores_bulkonly_top815.csv"

goi <- read.csv(paste(dir_goi, fn_goi, sep=""))

umi <- umi[rownames(umi) %in% goi$Symbol, ]
dim(umi)

# Order cells based on fine-level cell types (e.g., SC-F1)
meta <- arrange(meta, cluster)

umi <- umi[, rownames(meta)]

library(matrixStats)
row_vars<- rowVars(umi)

which(row_vars == 0)

rownames(umi)[which(row_vars==0)]

umi <- umi[-c(which(row_vars==0)), ]

#umi_hc <- pheatmap(umi, clustering_distance_rows="correlation", clustering_method="average", cluster_cols = F, scale = "row") 

df_celltypes <- data.frame(row.names = row.names(meta), Cell_Type = meta$cluster)

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- 'pain_genes_single_cell_pheatmap.png'

# umi_hc <- pheatmap(umi, color=colorRampPalette(c("purple", "black","yellow"))(50), clustering_distance_rows='correlation', show_colnames = F, show_rownames = F, cluster_cols = F, annotation_col = df_celltypes, main="Pain Associated Single Cell Heatmap in RA",  legend=T, legend_labels = c("row min zscore", "row max zscore"), file=paste(save_dir,fn,sep=""))
# umi_hc_rows <- rownames(umi[umi_hc$tree_row[["order"]],])

# TODO: we probably still need to use the DoHeatmap function. The expression matrix, even after row-wise rescaling, is way too sparse.
library(Seurat)
umi_seurat <- CreateSeuratObject(counts = umi, project = "normalized")
cell_types <- meta$cluster
umi_seurat@meta.data$groups <- as.factor(cell_types)

umi_seurat <- ScaleData(umi_seurat, features = rownames(umi_seurat))

umi_scaled <- umi_seurat@assays$RNA@scale.data

umi_hc <- pheatmap(umi_scaled, color=colorRampPalette(c("purple", "black","yellow"))(50), clustering_distance_rows='correlation', show_colnames = F, show_rownames = F, cluster_cols = F, annotation_col = df_celltypes, main="Pain Associated Single Cell Heatmap in RA",  legend=T, legend_labels = c("row min zscore", "row max zscore"), file=paste(save_dir,fn,sep=""))
#dev.off()
umi_hc_rows <- rownames(umi[umi_hc$tree_row[["order"]],])

# DoHeatmap(
#   umi_seurat,
#   #group.by = "groups",
#   features = umi_hc_rows
# )
# 
# dev.off()

heatmapfig <- DoHeatmap(
  umi_seurat,
  group.by = "groups",
  angle = 35,
  size = 2,
  features = umi_hc_rows
) + theme(axis.text.y = element_text(size = 0), axis.text.x = element_text(size = 0))

ggsave(filename = paste(save_dir, 'DoHeatmap_all_new_subtypes.pdf', sep=""), plot = heatmapfig)
#dev.copy(pdf,paste(save_dir, 'DoHeatmap_all_new_subtypes.pdf', sep=""))
#dev.off()







