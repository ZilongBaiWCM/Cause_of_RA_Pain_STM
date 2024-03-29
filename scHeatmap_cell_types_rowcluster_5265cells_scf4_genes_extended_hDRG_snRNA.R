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

m1 <- read.csv("/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/Ligand_Receptor/NicheNet_Results/RA_pain_ligands_to_hDRG_snRNA_receptors_interaction_interactome_all_neurons_vs_all_others.csv")
m1 <- m1[m1$ligand_type == 'SC-F4',]
scf4_ligands <- unique(m1$ligand)

m2 <- read.csv("/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/Ligand_Receptor/NicheNet_Results/hDRG_snRNA_ligands_to_RA_pain_receptors_interaction_interactome_all_neurons_vs_all_others.csv") 
m2 <- m2[m2$receptor_type == 'SC-F4',]
scf4_receptors <- unique(m2$receptor)

scf4_genes <- union(scf4_ligands, scf4_receptors)

umi <- umi[rownames(umi) %in% scf4_genes, ]
dim(umi)

# Order cells based on fine-level cell types (e.g., SC-F1)
meta <- arrange(meta, cluster)

umi <- umi[, rownames(meta)]

# Load DEA sorted gene list.
df_dea <- read.csv(paste('/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/', 'DEA_SC-F4_scf4_genes_sorted_genes_hDRG_snRNA.csv', sep=""), row.names = 1)
scf4_genes_dea <- rownames(df_dea)

# Rank rows of UMI according to DEA sorted gene list.
umi <- umi[scf4_genes_dea, ]

library(matrixStats)
row_vars<- rowVars(umi)

which(row_vars == 0)

rownames(umi)[which(row_vars==0)]

umi <- umi[which(row_vars!=0), ]
dim(umi)
#umi_hc <- pheatmap(umi, clustering_distance_rows="correlation", clustering_method="average", cluster_cols = F, scale = "row") 

df_celltypes <- data.frame(row.names = row.names(meta), Cell_Type = meta$cluster)

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- 'genes_scf4_single_cell_pheatmap.png'

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

print(dim(umi_seurat))

# DoHeatmap(
#   umi_seurat,
#   #group.by = "groups",
#   features = umi_hc_rows
# )
# 
# dev.off()

# heatmapfig <- DoHeatmap(
#   umi_seurat,
#   group.by = "groups",
#   angle = 35,
#   size = 2,
#   features = rownames(umi_seurat) #umi_hc_rows
# ) + theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 0))
# 
# ggsave(filename = paste(save_dir, 'DoHeatmap_all_new_subtypes_scf4_genes_DeaRanked_hDRG_snRNA.pdf', sep=""), plot = heatmapfig, dpi = 300 ) #, width = 2.2, height = 2.9, dpi = 180, units = "in")
# 
# width <- 2400  # Increased width
# height <- 1000  # Increased height
# 
# png(filename = paste(save_dir, 'DoHeatmap_all_new_subtypes_scf4_genes_DeaRanked_hDRG_snRNA.png', sep=""), width = width, height = height, res = 600)
# print(heatmapfig)
# dev.off()
#dev.copy(pdf,paste(save_dir, 'DoHeatmap_all_new_subtypes.pdf', sep=""))
#dev.off()







