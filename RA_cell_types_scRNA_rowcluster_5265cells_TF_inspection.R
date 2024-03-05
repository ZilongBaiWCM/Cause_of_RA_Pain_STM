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

tf_expression = umi['TF',]

meta$tf_expression <- tf_expression[rownames(meta)]

meta$cluster <- as.factor(meta$cluster)

# Define the color palette from RColorBrewer
#colors <- brewer.pal(n = length(unique(meta$cluster)), name = "Set1")
num_groups = length(unique(meta$cluster))
rainbow_colors <- rainbow(num_groups)

#brewer_rainbow_palette <- brewer.pal(n = num_groups, name = rainbow_colors)

# Create a boxplot with automatically set colors
TF_in_scRNAseq <- ggplot(meta, aes(x = cluster, y = tf_expression, fill = cluster)) +
  #geom_boxplot()+ #(outlier.shape = NA) +
  geom_jitter(aes(color = cluster), width = 0.2, size=3) +
  #coord_cartesian(ylim = c(0, 6)) +
  scale_fill_manual(values = rainbow_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text = element_text(size = 12)) +
  labs(x = "Cell Type", y = "log2(CPM+1) transformed UMIs counts")

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

ggsave(filename = paste(save_dir, 'transferrin_in_RA_synovium_scRNAseq.pdf', sep=""), plot = TF_in_scRNAseq, dpi = 300 ) #, width = 2.2, height = 2.9, dpi = 180, units = "in")


