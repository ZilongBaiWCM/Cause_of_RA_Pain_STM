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

data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"
gex_ENSG_Symbols <- read.csv(paste(data_dir, 'gex_ENSG_Symbols.csv', sep=""))
row.names(gex_ENSG_Symbols) <- gex_ENSG_Symbols$X

# Convert row names to a new column
sc_df <- sc_df %>% rownames_to_column()
gex_ENSG_Symbols <- gex_ENSG_Symbols %>% rownames_to_column()

# Perform left join based on row names
merged_df <- left_join(sc_df, gex_ENSG_Symbols, by = "rowname")

tf_expression = merged_df[complete.cases(merged_df) & merged_df$Symbol=='TF',!(names(merged_df) %in% c('rowname', 'X', 'ENSG', 'Symbol'))]


# Load meta information for the sorted cell expression column (i.e., tube) groups.
file = paste(sc_dir,"low_input_meta.725715.tsv", sep="")
sc_meta <- read.delim2(file, header = TRUE, sep = "\t", dec = ".")

rownames(sc_meta) <- sc_meta$Tube.label

common_tubes <- intersect(names(tf_expression), row.names(sc_meta))

sc_meta <- sc_meta[common_tubes,]
tf_expression <- tf_expression[,common_tubes]

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'



#-------------------------- QC..pass.fail. and Tissue.type RA --------------------------#

attach(sc_meta)
sc_meta_qc <- sc_meta[which(QC..pass.fail.=="passed"), ]
detach(sc_meta)

attach(sc_meta_qc)
meta <- sc_meta_qc[which(Tissue.type %in% c("RA-biopsy", "RA-arthro")),]
detach(sc_meta_qc)

tf_expression <- tf_expression[,rownames(meta)]

meta$tf_expression <- t(tf_expression[rownames(meta)])

meta$Cell.type <- as.factor(meta$Cell.type)

num_groups = length(unique(meta$Cell.type))
rainbow_colors <- rainbow(num_groups)

# Create a boxplot with automatically set colors
TF_in_sorted <- ggplot(meta, aes(x = Cell.type, y = tf_expression, fill = Cell.type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Cell.type), width = 0.2, size=1.3) +
  #coord_cartesian(ylim = c(0, 6)) +
  scale_fill_manual(values = rainbow_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text = element_text(size = 7),
        legend.position = "none") +
  labs(x = "Cell Type", y = "Transcripts Per Million (TPMs)")

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

library(svglite)
#ggsave(filename = paste(save_dir, 'transferrin_in_RA_synovium_sorted_bulk_RNAseq.pdf', sep=""), plot = TF_in_sorted, dpi = 300 ) #, width = 2.2, height = 2.9, dpi = 180, units = "in")
ggsave(filename = paste(save_dir, 'transferrin_in_RA_synovium_sorted_bulk_RNAseq.svg', sep=""), plot = TF_in_sorted, device="svg", width = 3.6, height = 3)

