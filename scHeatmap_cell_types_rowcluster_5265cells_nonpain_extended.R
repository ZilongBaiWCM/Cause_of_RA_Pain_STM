#' ---
#' title: "Single-cell RNA-seq clustering"
#' Integrative pipeline of integrating bulk RNA-seq with single-cell RNA-seq and then unbiased clustering
#' author: "Fan Zhang"
#' date: "2018-03-19"
#' ---
#' 
#setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")
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

# -------
# Main steps
# 1) We first selected the highly variable genes from both sides that we 
#    can keep the data set specific variation; 
# 2) Based on the shared highly variable genes, we integrate single-cell RNA-seq 
#    with bulk RNA-seq by CCA (González et al. 2008); 
# 3) We calculate the cell-to-cell similarity matrix on top 10 CCA canonical 
#    dimensions based on Euclidean distance;
# 4) We build up a K-nearest neighbors (KNN) graph based on the cell-to-cell similarity matrix; 
#    We then convert the KNN neighbor relation matrix into an adjacency matrix; 
# 5) We cluster the cells using Infomap algorithm (cluster_infomap function from igraph package) 
#    for community detection to decompose the cell-to-cell adjacency matrix into major modules; 
# 6) We then construct a low dimensional embedding using tSNE based on 
#    the cell-to-cell distance matrix with cell clusters labeled; 
# 7) We identify and prioritize significantly differential expressed (DE) genes for 
#    each distinct cluster using AUC and Wilcox test; 
# 8) For pathway analysis, we downloaded gene sets from Gene Ontology (GO) terms on April 2016. 
#    We also use the immunologic signatures from 4,872 hallmark gene sets from MSigDB (Liberzon et al. 2015)
#    to test enrichment of all the tested DE genes sorted by decreased AUC scores 
#    for each cluster using liger.
#
# NOTE: 
# 1) coarse clustering and fine clustering follow the same pipeline
#    For fine clustering, please run the piepeline on cells from one cell type 
# 2) Since graph-based clustering is based on random walk (a stochastic process), the results depand on nk (KNN step), 
#    and also the number dimensions ncc (CCA step), and the number of selected features as well.
# 
# -------


# Read the post-QC single-cell RNA-seq expression data and meta data
#dat <- readRDS(file = paste("../data/celseq_synovium_log2_postQC", ".rds", sep = ""))
## NOTE: We are now using the normalized data provided by Dr. Fan Zhang. 
## We compute the normalized counts with different normalization functions from the raw counts
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
# Pain_associated genes
dir_pain<- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/"
fn_pain <- "Laplacian_Scores_of_zscores_bulkonly_top815.csv"
genes_pain <- read.csv(paste(dir_pain, fn_pain, sep=""))

# Low-inflammatory genes
data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"
gex_ENSG_Symbols <- read.csv(paste(data_dir, 'gex_ENSG_Symbols.csv', sep=""))
row.names(gex_ENSG_Symbols) <- gex_ENSG_Symbols$X

ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores', filter_condition,'.csv', sep=""), row.names=1)

genes_lowinf <- gex_ENSG_Symbols[rownames(ls_result), 'Symbol']

genes_nonpain <- genes_lowinf[!(genes_lowinf %in% genes_pain$Symbol)]

umi <- umi[rownames(umi) %in% genes_nonpain, ]
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

fn <- 'nonpain_genes_single_cell_pheatmap.png'

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

ggsave(filename = paste(save_dir, 'DoHeatmap_all_new_subtypes_nonpain.pdf', sep=""), plot = heatmapfig)
#dev.copy(pdf,paste(save_dir, 'DoHeatmap_all_new_subtypes.pdf', sep=""))
#dev.off()

# TODO: 
# 1. align the columns to specific cell types. Is it done? Yes, actually. Need to put them on top of the columns to indicate.
# 2. need to access the visualization problem. Why are they not visualized in a correct way?
# some rows are zeros?
# also, why is the data not scaled?
# dev.off()





