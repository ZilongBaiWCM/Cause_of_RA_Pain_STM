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
# sc_dir <- "/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/Fibroblast_subtyping/amp_phase1_ra-master/data/"
# fn_nd = paste(sc_dir,"celseq_synovium_log2_5265cells_paper.rds", sep="")
# 
# fn_meta = paste(sc_dir,"celseq_synovium_meta_5265cells_paper.rds", sep="")
# 
# umi_5265cells <- readRDS(fn_nd)
#   
# meta_5265cells <- readRDS(fn_meta)
# rownames(meta_5265cells) <- meta_5265cells$cell_name
# 
# meta_RA <- meta_5265cells[meta_5265cells$disease=="RA", ]
# umi_RA <- umi_5265cells[, rownames(meta_RA)]
# 
# umi <- umi_RA # Only focus on RA cells in this analysis.
# meta <- meta_RA

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

#write.table(genes_nonpain,file='/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/nonpain_low_inflammation_gene_symbols.csv',col.names=FALSE, row.names=FALSE)

#write.table(gex_ENSG_Symbols$Symbol,file='/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/bulk_background_gene_symbols.csv',col.names=FALSE, row.names=FALSE)

painVSnonpain <- gex_ENSG_Symbols[rownames(ls_result), ]
painVSnonpain$pain_status <- 'Nonpain'
painVSnonpain[painVSnonpain$Symbol %in% genes_pain$Symbol, 'pain_status'] <- 'Pain'

write.csv(painVSnonpain,file='/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/low_inflammatory_painVSnonpain_genes.csv')
