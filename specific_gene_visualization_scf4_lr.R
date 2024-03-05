#' ---
#' title: 
#' 
#' author: "Zilong Bai"
#' date: "2023-02-23"
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
library(biomaRt)
#library(edgeR) # To use the cpm function for computing log2cpm from integer counts.
require(gdata)

source("meta_colors.R")
source("pure_functions.R")




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

log2cpm <- umi_RA # Only focus on RA cells in this analysis.
meta <- meta_RA

all(meta$cell_name == colnames(log2cpm))
log2cpm <- as.data.frame(log2cpm)


m1 <- read.csv("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/final_hdrg_tables_and_code/RA_ligand_hDRG_receptor_pairs_Seurat_Fibroblast_vsOthers_padj_genes_f4_yes_c19.csv")
scf4_ligands <- unique(m1$ligand_gname)

m2 <- read.csv("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/final_hdrg_tables_and_code/hDRG_ligand_RA_receptor_pairs_Seurat_Fibroblast_vsOthers_padj_genes_f4_yes_c19.csv") 
scf4_receptors <- unique(m2$receptor_gname)

scf4_genes <- union(scf4_ligands, scf4_receptors)


fontsize <- 30

# -----
# Setting color palattee
mycolors = c(brewer.pal(name="RdBu", n = 8), brewer.pal(name="Oranges", n = 4), brewer.pal(name="BrBG", n = 6))

# -----
#save_dir <- '/Users/baizilong/Documents/Dana-dataset/Single_Cell_Analysis/'
save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

for(marker in scf4_genes){
  meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == marker),])
  dat_percent <- meta %>%
    group_by(cluster) %>%
    summarise(percent = sum(marker > 0) / length(marker) * 100)
  
  p1 <- ggplot() +
    geom_col(
      data = dat_percent,
      mapping = aes(x = forcats::fct_rev(cluster), y = percent, fill = cluster),
      color = "grey50", 
      size = 0.1
    ) +
    scale_fill_manual(values=mycolors)+
    scale_x_discrete(limits = rev(levels(meta$cluster))) +
    #scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
    #scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    ylim(0,100)+
    geom_vline(xintercept = 14.5, linetype="dashed") +
    geom_vline(xintercept = 10.5, linetype="dashed") +
    geom_vline(xintercept = 6.5, linetype="dashed") +
    coord_flip() +
    labs(x = NULL, y = "% non-zero", title = marker) +
    theme_bw(base_size = fontsize) + theme(
      legend.position = "none",
      #axis.text.y       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title = element_text(color="black", size=fontsize, face = "italic")
    )
  
  # data_exp <- data.frame(matrix(ncol = 0, nrow = sum(meta$marker>0)))
  # data_exp$gex <- meta$marker[meta$marker>0]
  # data_exp$cluster <- meta$cluster[meta$marker>0]
  data_exp <- data.frame(matrix(ncol = 0, nrow = length(meta$marker)))
  data_exp$gex <- meta$marker
  data_exp$cluster <- meta$cluster
  p2 <- ggplot(data_exp, aes(x = forcats::fct_rev(cluster), y = gex, color = cluster)) +
    geom_jitter()+
    # geom_violin(
    #   data = data_exp, #dat_percent,
    #   mapping = aes(x = forcats::fct_rev(cluster), y = gex, fill = cluster),
    #   color = "grey50", size = 0.2
    # ) +
    # geom_jitter(data = data_exp, #dat_percent,
    #             mapping = aes(x = forcats::fct_rev(cluster), y = gex, fill = cluster),
    #             #color = "grey50", 
    #             size = 0.3,
    #             position = position_jitter(seed = 1, width = 0.2)) +
    scale_x_discrete(limits = rev(levels(meta$cluster))) +
    scale_color_manual(values=mycolors) +
    #scale_color_brewer(palette="RdBu") + 
    #scale_fill_manual(values = meta_colors$cluster, name = "Cluster") +
    #scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    geom_vline(xintercept = 14.5, linetype="dashed") +
    geom_vline(xintercept = 10.5, linetype="dashed") +
    geom_vline(xintercept = 6.5, linetype="dashed") +
    coord_flip() +
    labs(x = NULL, y = expression("Log"['2']*"(CPM+1)"), title = NULL) +
    theme_bw(base_size = fontsize) + theme(
      legend.position = "none",
      axis.text.y       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      panel.border    = element_rect(size = 0.5),
      plot.title = element_text(color="black", size=fontsize, face = "italic")
    )
  
  
  g_ntn4 <- egg::ggarrange(
    p1, p2,
    ncol = 2
  )
  
  ggsave(paste(marker,'_figure6A_scf4_genes.png',sep=""), plot=g_ntn4, path=save_dir, width = 8, height = 10, dpi = 180, units = "in", device='png')
  rm(p1)
  rm(p2)
}


# ============================================== Fin =================================================== #



