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
