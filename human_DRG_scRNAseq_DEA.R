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
library(biomaRt)
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

m1 <- read.csv("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/final_hdrg_tables_and_code/RA_ligand_hDRG_receptor_pairs_Seurat_Fibroblast_vsOthers_padj_genes_f4_yes_c19.csv")
scf4_ligands <- unique(m1$ligand_gname)

m2 <- read.csv("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/final_hdrg_tables_and_code/hDRG_ligand_RA_receptor_pairs_Seurat_Fibroblast_vsOthers_padj_genes_f4_yes_c19.csv") 
scf4_receptors <- unique(m2$receptor_gname)

scf4_genes <- union(scf4_ligands, scf4_receptors)

umi <- umi[rownames(umi) %in% scf4_genes, ]
dim(umi)

# Only focus on SC-F cells:
meta <- meta[meta$cluster %in% c('SC-F1', 'SC-F2', 'SC-F3', 'SC-F4'), ]

# Order cells based on fine-level cell types (e.g., SC-F1)
meta <- arrange(meta, cluster)

umi <- umi[, rownames(meta)]

library(matrixStats)
row_vars<- rowVars(umi)

which(row_vars == 0)

rownames(umi)[which(row_vars==0)]

umi <- umi[which(row_vars!=0), ]

#umi_hc <- pheatmap(umi, clustering_distance_rows="correlation", clustering_method="average", cluster_cols = F, scale = "row") 

df_celltypes <- data.frame(row.names = row.names(meta), Cell_Type = meta$cluster)

#save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'
#'/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Plots_for_paper/'


df_celltypes <- df_celltypes[colnames(umi), 1, drop=FALSE]

colnames(umi) <- paste(df_celltypes$Cell_Type, colnames(umi),sep="_")

library(Seurat)
umi_seurat <- CreateSeuratObject(counts = umi, project = "normalized")
cell_types <- meta$cluster
umi_seurat@meta.data$groups <- as.factor(cell_types)

# NOTE: This seurat object is built with counts that have already been log normalized. The later on operations on the data slot is on the normalized data.

cell_types <- levels(umi_seurat)
de.markers.list <- list()
ctype.de.markers.list <- list()

library(ggrepel)

