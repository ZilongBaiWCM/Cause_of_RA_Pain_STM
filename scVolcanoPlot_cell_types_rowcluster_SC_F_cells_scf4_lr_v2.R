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

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'
#'/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Plots_for_paper/'

# umi_hc <- pheatmap(umi, color=colorRampPalette(c("purple", "black","yellow"))(50), clustering_distance_rows='correlation', show_colnames = F, show_rownames = F, cluster_cols = F, annotation_col = df_celltypes, main="Pain Associated Single Cell Heatmap in RA",  legend=T, legend_labels = c("row min zscore", "row max zscore"), file=paste(save_dir,fn,sep=""))
# umi_hc_rows <- rownames(umi[umi_hc$tree_row[["order"]],])


df_celltypes <- df_celltypes[colnames(umi), 1, drop=FALSE]

colnames(umi) <- paste(df_celltypes$Cell_Type, colnames(umi),sep="_")

library(Seurat)
umi_seurat <- CreateSeuratObject(counts = umi, project = "normalized")
cell_types <- meta$cluster
umi_seurat@meta.data$groups <- as.factor(cell_types)

# NOTE: This seurat object is built with counts that have already been log normalized. The later on operations on the data slot is on the normalized data.

cell_types_levels <- levels(umi_seurat)
de.markers.list <- list()
ctype.de.markers.list <- list()

library(ggrepel)

for (ctype in cell_types_levels) {
  de.markers <- FindMarkers(umi_seurat, ident.1 = ctype, ident.2 = NULL, logfc.threshold = 0,
                            min.pct = 0,
                            min.cells.feature = 0,
                            min.cells.group = 0,
  )
  print(dim(de.markers))
  ctype.de.markers = de.markers[de.markers$avg_log2FC > 0 & de.markers$p_val_adj < 0.05, ]
  print(dim(ctype.de.markers))
  
  de.markers.list[[length(de.markers.list)+1]] <- de.markers
  ctype.de.markers.list[[length(ctype.de.markers.list)+1]] <- ctype.de.markers
  
  # add a column of NAs
  de.markers[['DiffExp']] <- "NO"
  # if log2Foldchange > 0 and pvalue < 0.05, set as "UP" 
  de.markers$DiffExp[de.markers$avg_log2FC > 0 & de.markers$p_val_adj < 0.05] <- "UP"
  # if log2Foldchange < 0 and pvalue < 0.05, set as "DOWN"
  de.markers$DiffExp[de.markers$avg_log2FC < 0 & de.markers$p_val_adj < 0.05] <- "DOWN"
  
  color_vals <- c()
  if("DOWN" %in% de.markers$DiffExp){
    color_vals <- append(color_vals, "blue")
  }
  if("NO" %in% de.markers$DiffExp){
    color_vals <- append(color_vals, "grey")
  }
  if("UP" %in% de.markers$DiffExp){
    color_vals <- append(color_vals, "red")
  }
  
  g1 <- de.markers[rownames(de.markers) == "NTN4", ]
  g2 <- de.markers[rownames(de.markers) == "HBEGF", ]
  
  p<-ggplot(data=de.markers, aes(x=avg_log2FC, y=-log(p_val_adj), col=DiffExp)) +
    geom_point(aes(size=1.2), show.legend=FALSE, size=1.2) + 
    theme(text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size=12)) +
    theme_minimal() +
    scale_color_manual(values=color_vals) +
    geom_vline(xintercept=c(0), col="red") +
    geom_hline(yintercept=-log(0.05), col="red") +
    ylim(-10,450) +
    xlim(-10, 10) +
    labs(x=NULL, y=NULL)
  p <- p + geom_point(aes(size=1.2), data=g1, shape=1, show.legend=FALSE) +
    geom_text(data=g1, label="NTN4", vjust=-1, size=3) +
    labs(x=NULL, y=NULL) +
    theme(legend.position = "none")
  p + geom_point(aes(size=1.2), data=g2, shape=1, show.legend=FALSE) +
    geom_text(data=g2, label="HBEGF", vjust=-1, size=3) +
    labs(x=NULL, y=NULL) +
    theme(legend.position = "none")
    
  
  fn <- paste(save_dir, 'DEA_',ctype,'_volcano_plot_figure3_scf4_genes_extended.png',sep="")
  ggsave(file=fn, bg='white', width = 2.2, height = 2.9, dpi = 180, units = "in", device='png')
  
  sorted_genes <- de.markers[order(-de.markers$avg_log2FC), ]
  write.csv(sorted_genes, paste(save_dir, 'DEA_',ctype,'_scf4_genes_sorted_genes.csv',sep=""))
  #dev.off()
}
