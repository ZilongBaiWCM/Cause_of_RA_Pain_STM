# This script aims to find pathways that are "similar" to the Neuroactive Ligand-Receptor Interaction.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
rm(list=ls())
library(magrittr)
library(dplyr)
library(kim)
library(pheatmap)
library(RColorBrewer)
#BiocManager::install("M3C")

#library(M3C)

#help(clustersim)
#help(M3C)

input_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Generated_Gene_lists/"

df_genes <- read.csv(paste(input_dir, 'bulk_padj_<0.01log2FC_<0.csv',sep=""))
genes_lowinf <- df_genes$X
rm(df_genes)

# Load genes collectively significantly correlate with pain scores in bulk RNA-seq data.
ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
s <- 815
file<-'Laplacian_Scores_of_zscores_bulkonly_top815_in_sorted_769.csv'
laplacian_score_top <- read.csv(paste(ls_dir,file,sep=""))
row.names(laplacian_score_top) <- laplacian_score_top$X

laplacian_score_top <- laplacian_score_top[,2:4]

# Load sorted cell expression data.
sc_dir <- "/Users/baizilong/Documents/Dana-dataset/Bulk Sorted/"

file = paste(sc_dir,"low_input_gene_sample_tpm_matrix.725714.tsv", sep="")
sc_df <- read.delim2(file, header = TRUE, sep = "\t", dec = ".")
dim(sc_df)
rownames(sc_df) <- sc_df$X
sc_df <- sc_df[,-1]

# Extract sorted cell expression of the low inflammatory genes that are not pain-associated.
genes_nonpain <- genes_lowinf[!(genes_lowinf %in% row.names(laplacian_score_top))]
sc_nonpain_in_bulk <- sc_df[row.names(sc_df) %in% genes_nonpain,]

# Load meta information for the sorted cell expression column (i.e., tube) groups.
file = paste(sc_dir,"low_input_meta.725715.tsv", sep="")
sc_meta <- read.delim2(file, header = TRUE, sep = "\t", dec = ".")

rownames(sc_meta) <- sc_meta$Tube.label

common_tubes <- intersect(colnames(sc_nonpain_in_bulk), row.names(sc_meta))

sc_meta <- sc_meta[common_tubes,]
sc_nonpain_in_bulk <- sc_nonpain_in_bulk[,common_tubes]

# Regroup the tubes based on their corresponding cell types.
sc_meta_tube_celltypes <- data.frame(row.names = row.names(sc_meta), Tubes = row.names(sc_meta), Cell_Type = sc_meta$Cell.type)

attach(sc_meta_tube_celltypes)
sc_meta_tube_celltypes <- sc_meta_tube_celltypes[order(Cell_Type),]
detach(sc_meta_tube_celltypes)

df_celltypes <- data.frame(row.names = row.names(sc_meta_tube_celltypes), Cell_Type = sc_meta_tube_celltypes$Cell_Type)

sc_nonpain_in_bulk <- sc_nonpain_in_bulk[, row.names(df_celltypes)]

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- "Sorted Cell Expressions of Genes Collectively Correlate with Pain Scores nonpain_in_sorted.png"
#pheatmap(sc_nonpain_in_bulk, color=colorRampPalette(c("blue", "white","red"))(50), scale="row", show_colnames = F, show_rownames = F, cluster_cols = F, annotation_col = df_celltypes, main="Sorted Cell Expressions of Genes Collectively Correlate with Pain Scores",  legend=T, legend_labels = c("row min zscore", "row max zscore"), file=paste(save_dir,fn,sep=""))



#-------------------------- QC..pass.fail. and Tissue.type RA --------------------------#

attach(sc_meta)
sc_meta_qc <- sc_meta[which(QC..pass.fail.=="passed"), ]
detach(sc_meta)

attach(sc_meta_qc)
sc_meta_qc_tt <- sc_meta_qc[which(Tissue.type %in% c("RA-biopsy", "RA-arthro")),]
detach(sc_meta_qc)

sc_nonpain_in_bulk_qc_tt <- sc_nonpain_in_bulk[,rownames(sc_meta_qc_tt)]

# Regroup the tubes based on their corresponding cell types.
sc_meta_tube_celltypes_qc_tt <- data.frame(row.names = row.names(sc_meta_qc_tt), Tubes = row.names(sc_meta_qc_tt), Cell_Type = sc_meta_qc_tt$Cell.type)

attach(sc_meta_tube_celltypes_qc_tt)
sc_meta_tube_celltypes_qc_tt <- sc_meta_tube_celltypes_qc_tt[order(Cell_Type),]
detach(sc_meta_tube_celltypes_qc_tt)

df_celltypes_qc_tt <- data.frame(row.names = row.names(sc_meta_tube_celltypes_qc_tt), Cell_Type = sc_meta_tube_celltypes_qc_tt$Cell_Type)

df_celltypes_qc_tt <- df_celltypes_qc_tt %>% rename(`Cell Type` = Cell_Type)

sc_nonpain_in_bulk_qc_tt <- sc_nonpain_in_bulk_qc_tt[, row.names(df_celltypes_qc_tt)]

ra_passed_sc_nonpain_in_bulk_qc_tt <- sc_nonpain_in_bulk_qc_tt

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- "Bulk Sorted Cell QC passed RA Tissue nonpain_in_sorted_1458.png"
# pheatmap(sc_nonpain_in_bulk_qc_tt, 
#          color=colorRampPalette(c("green", "white","red"))(50), 
#          scale="row", 
#          show_colnames = F, 
#          show_rownames = F, 
#          cluster_cols = F, 
#          treeheight_row = 0,
#          annotation_col = df_celltypes_qc_tt, 
#          fontsize = 12,
#          #main="Bulk Sorted Cell Expressions (QC passed, RA tissue)",  
#          legend=T, legend_labels = c("row min zscore", "row max zscore"), 
#          file=paste(save_dir,fn,sep=""))

pheatmap(sc_nonpain_in_bulk_qc_tt, 
         color=colorRampPalette(c("blue", "white","red"))(50), 
         scale="row", 
         show_colnames = F, 
         show_rownames = F, 
         cluster_cols = F, 
         treeheight_row = 0,
         cellheight=.5, cellwidth = 6,
         annotation_col = df_celltypes_qc_tt, 
         fontsize = 12,
         #main="Bulk Sorted Cell Expressions (QC passed, RA tissue)",  
         legend=T, legend_labels = c("row min zscore", "row max zscore"),
         file=paste(save_dir,fn,sep=""))


fn <- "GEX Bulk Sorted Cell QC passed RA Tissue nonpain_in_sorted.csv"
write.csv(sc_nonpain_in_bulk_qc_tt, paste(save_dir,fn,sep=""))

