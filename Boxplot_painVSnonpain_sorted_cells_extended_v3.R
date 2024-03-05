# This script aims to find pathways that are "similar" to the Neuroactive Ligand-Receptor Interaction.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
rm(list=ls())
library(magrittr)
library(dplyr)
library(kim)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
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

# Extract pain-associated genes in sorted cell GEX
genes_pain <- row.names(laplacian_score_top)

# Load sorted cell expression data.
sc_dir <- "/Users/baizilong/Documents/Dana-dataset/Bulk Sorted/"

file = paste(sc_dir,"low_input_gene_sample_tpm_matrix.725714.tsv", sep="")
sc_df <- read.delim2(file, header = TRUE, sep = "\t", dec = ".")
dim(sc_df)
rownames(sc_df) <- sc_df$X
sc_df <- sc_df[,-1]

# Extract non-pain-associated low-inflammatory genes in sorted cell GEX
genes_nonpain <- genes_lowinf[!(genes_lowinf %in% genes_pain)]

# Extract sorted cell expression of the pain-associated genes.
sc_pain_in_bulk <- sc_df[row.names(sc_df) %in% genes_pain,]
# Extract sorted cell expression of the low inflammatory genes that are not pain-associated.
sc_nonpain_in_bulk <- sc_df[row.names(sc_df) %in% genes_nonpain,]

# Load meta information for the sorted cell expression column (i.e., tube) groups.
file = paste(sc_dir,"low_input_meta.725715.tsv", sep="")
sc_meta <- read.delim2(file, header = TRUE, sep = "\t", dec = ".")

rownames(sc_meta) <- sc_meta$Tube.label

common_tubes <- intersect(colnames(sc_nonpain_in_bulk), row.names(sc_meta))

sc_meta <- sc_meta[common_tubes,]

sc_pain_in_bulk <- sc_pain_in_bulk[,common_tubes]
sc_nonpain_in_bulk <- sc_nonpain_in_bulk[,common_tubes]

# Regroup the tubes based on their corresponding cell types.
sc_meta_tube_celltypes <- data.frame(row.names = row.names(sc_meta), Tubes = row.names(sc_meta), Cell_Type = sc_meta$Cell.type)

attach(sc_meta_tube_celltypes)
sc_meta_tube_celltypes <- sc_meta_tube_celltypes[order(Cell_Type),]
detach(sc_meta_tube_celltypes)

df_celltypes <- data.frame(row.names = row.names(sc_meta_tube_celltypes), Cell_Type = sc_meta_tube_celltypes$Cell_Type)

sc_pain_in_bulk <- sc_pain_in_bulk[, row.names(df_celltypes)]
sc_nonpain_in_bulk <- sc_nonpain_in_bulk[, row.names(df_celltypes)]


#-------------------------- QC..pass.fail. and Tissue.type RA --------------------------#

attach(sc_meta)
sc_meta_qc <- sc_meta[which(QC..pass.fail.=="passed"), ]
detach(sc_meta)

attach(sc_meta_qc)
sc_meta_qc_tt <- sc_meta_qc[which(Tissue.type %in% c("RA-biopsy", "RA-arthro")),]
detach(sc_meta_qc)

attach(sc_meta_qc_tt)
sc_meta_qc_tt <- sc_meta_qc_tt[order(Cell.type),]
detach(sc_meta_qc_tt)

sc_pain_in_bulk_qc_tt <- sc_pain_in_bulk[,sc_meta_qc_tt$Tube.label]
sc_nonpain_in_bulk_qc_tt <- sc_nonpain_in_bulk[,sc_meta_qc_tt$Tube.label]

df_celltypes<- data.frame(row.names = sc_meta_qc_tt$Tube.label, Tube_Label = sc_meta_qc_tt$Tube.label,Cell_Type = sc_meta_qc_tt$Cell.type)


attach(df_celltypes)
df_celltypes <- df_celltypes[order(Cell_Type),]
detach(df_celltypes)

sc_pain_in_bulk_qc_tt <- sc_pain_in_bulk[,rownames(df_celltypes)]
sc_nonpain_in_bulk_qc_tt <- sc_nonpain_in_bulk[,rownames(df_celltypes)]

#--------------------------- Z-scoring gene expression --------------------------#
tpm_pain <- sc_pain_in_bulk_qc_tt

tpm_nonpain <- sc_nonpain_in_bulk_qc_tt

#--------------------------- Create data frame from GEX z-scores for boxplots --------------------------#
unique_cell_types <- unique(df_celltypes$Cell_Type)

cell_types_vec <- c()
pnp_vec <- c()
mean_zscore_vec <- c()
gene_names <- c()

for(ct in unique_cell_types){
  tubes <- rownames(df_celltypes[df_celltypes$Cell_Type==ct, ])
  pain_gex <- tpm_pain[,tubes]
  nonpain_gex <- tpm_nonpain[,tubes]
  
  mean_pain_gex <- rowMeans(pain_gex)
  mean_nonpain_gex <- rowMeans(nonpain_gex)
  
  mean_zscore_vec <- append(mean_zscore_vec, unname(mean_nonpain_gex))
  mean_zscore_vec <- append(mean_zscore_vec, unname(mean_pain_gex))
  
  pnp_vec <- append(pnp_vec, rep('Non-pain', length(mean_nonpain_gex)))
  pnp_vec <- append(pnp_vec, rep('Pain', length(mean_pain_gex)))
  
  gene_names <- append(gene_names, names(mean_nonpain_gex))
  gene_names <- append(gene_names, names(mean_pain_gex))
  
  cell_types_vec <- append(cell_types_vec, rep(ct, length(mean_nonpain_gex) + length(mean_pain_gex)))
  
  
}

data=data.frame(gene_names, cell_types_vec, pnp_vec ,  mean_zscore_vec)

#--------------------------- Generate boxplots and save --------------------------#

# grouped boxplot
p <- ggplot(data, aes(x=cell_types_vec, y=mean_zscore_vec, fill=pnp_vec)) + 
  geom_boxplot() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=22),
  ) +
  scale_fill_manual(values = c("cyan", "orange")) +
  labs(fill = "")

save_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

fn <- "Bulk Sorted Cell QC passed RA Tissue pain vs nonpain boxplots v2_v3.png"

# ggsave(plot = p, 
#        filename=paste(save_dir, fn,vsep=""),
#        width = 12, height = 6, dpi = 180, units = "in", device='png')

#--------------------------- T test compare pain vs non-pain per cell type --------------------------#

results<-c()
for(ct in unique_cell_types){
  data_ct <- data[data$cell_types_vec==ct, ]
  data_pain <- data_ct[data_ct$pnp_vec=='Pain', "mean_zscore_vec"]
  data_nonpain <- data_ct[data_ct$pnp_vec=='Non-pain', "mean_zscore_vec"]
  test_result <- t.test(data_pain, data_nonpain)
  pval<-test_result$p.value
  print(paste(ct, "Pain vs Non-pain T test p-value:", pval))
  results <- append(results, -log(pval))
}
results[is.infinite(results)] <- NA
names(results) <- unique_cell_types

data_forplot <- data.frame(Category = names(results), Value = results)

p<- ggplot(data_forplot, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "Pain vs Non-pain T test", x = "Cell Type", y = "-log p-value")

ggsave(plot = p, 
       filename=paste(save_dir, "Pain vs Non-pain T test_v3.png",sep=""),
       width = 12, height = 6, dpi = 180, units = "in", device='png')

results<-c()
for(ct in unique_cell_types){
  
  # Compute mean z-score of pain in ct
  tubes <- rownames(df_celltypes[df_celltypes$Cell_Type==ct, ])
  mean_pain_gex_ct <- rowMeans(tpm_pain[, tubes])
  nontubes <- rownames(df_celltypes[!(df_celltypes$Cell_Type==ct), ])
  mean_pain_gex_nonct <- rowMeans(tpm_pain[, nontubes])
  test_result_ctVSothers <- t.test(mean_pain_gex_ct, mean_pain_gex_nonct, paired = T)
  pval<-test_result_ctVSothers$p.value
  print(paste(ct, "vs Others Pain Paired T test p-value:", pval))
  results <- append(results, -log(pval))
}
results[is.infinite(results)] <- NA
names(results) <- unique_cell_types

data_forplot <- data.frame(Category = names(results), Value = results)

p<- ggplot(data_forplot, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "One vs Others Pain Paired T test", x = "Cell Type", y = "-log p-value")

ggsave(plot = p, 
       filename=paste(save_dir, "One vs Others Pain Paired T test_v3.png",sep=""),
       width = 12, height = 6, dpi = 180, units = "in", device='png')

results<-c()
for(ct in unique_cell_types){
  
  # Compute mean z-score of pain in ct
  tubes <- rownames(df_celltypes[df_celltypes$Cell_Type==ct, ])
  mean_nonpain_gex_ct <- rowMeans(tpm_nonpain[, tubes])
  nontubes <- rownames(df_celltypes[!(df_celltypes$Cell_Type==ct), ])
  mean_nonpain_gex_nonct <- rowMeans(tpm_nonpain[, nontubes])
  test_result_ctVSothers <- t.test(mean_nonpain_gex_ct, mean_nonpain_gex_nonct, paired = T)
  pval<-test_result_ctVSothers$p.value
  print(paste(ct, "vs Others Non-Pain Paired T test p-value:", pval))
  results <- append(results, -log(pval))
}
results[is.infinite(results)] <- NA
names(results) <- unique_cell_types

data_forplot <- data.frame(Category = names(results), Value = results)

p<- ggplot(data_forplot, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "One vs Others Non-Pain Paired T test", x = "Cell Type", y = "-log p-value")

ggsave(plot = p, 
       filename=paste(save_dir, "One vs Others Non-Pain Paired T test_v3.png",sep=""),
       width = 12, height = 6, dpi = 180, units = "in", device='png')

#--------------------------- T test compare pain vs non-pain per cell type --------------------------#
results<-c()
for(ct in unique_cell_types){
  data_ct <- data[data$cell_types_vec==ct, ]
  data_pain <- data_ct[data_ct$pnp_vec=='Pain', "mean_zscore_vec"]
  data_nonpain <- data_ct[data_ct$pnp_vec=='Non-pain', "mean_zscore_vec"]
  test_result <- t.test(data_pain, data_nonpain, alternative="greater")
  pval<-test_result$p.value
  print(paste(ct, "Pain vs Non-pain T test Alternative greater p-value:", pval))
  results <- append(results, -log(pval))
}
results[is.infinite(results)] <- NA
names(results) <- unique_cell_types

data_forplot <- data.frame(Category = names(results), Value = results)

p<- ggplot(data_forplot, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "Pain vs Non-pain T test Alternative greater", x = "Cell Type", y = "-log p-value")

ggsave(plot = p, 
       filename=paste(save_dir, "Pain vs Non-pain T test Alternative greater_v3.png",sep=""),
       width = 12, height = 6, dpi = 180, units = "in", device='png')

results<-c()
for(ct in unique_cell_types){
  
  # Compute mean z-score of pain in ct
  tubes <- rownames(df_celltypes[df_celltypes$Cell_Type==ct, ])
  mean_pain_gex_ct <- rowMeans(tpm_pain[, tubes])
  nontubes <- rownames(df_celltypes[!(df_celltypes$Cell_Type==ct), ])
  mean_pain_gex_nonct <- rowMeans(tpm_pain[, nontubes])
  test_result_ctVSothers <- t.test(mean_pain_gex_ct, mean_pain_gex_nonct, paired = T, alternative="greater")
  pval<-test_result_ctVSothers$p.value
  print(paste(ct, "vs Others Pain Paired T test Alternative greater p-value:", pval))
  results <- append(results, -log(pval))
}
results[is.infinite(results)] <- NA
names(results) <- unique_cell_types

data_forplot <- data.frame(Category = names(results), Value = results)

p<- ggplot(data_forplot, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "One vs Others Pain Paired T test Alternative greater", x = "Cell Type", y = "-log p-value")

ggsave(plot = p, 
       filename=paste(save_dir, "One vs Others Pain Paired T test Alternative greater_v3.png",sep=""),
       width = 12, height = 6, dpi = 180, units = "in", device='png')

results<-c()
for(ct in unique_cell_types){
  
  # Compute mean z-score of pain in ct
  tubes <- rownames(df_celltypes[df_celltypes$Cell_Type==ct, ])
  mean_nonpain_gex_ct <- rowMeans(tpm_nonpain[, tubes])
  nontubes <- rownames(df_celltypes[!(df_celltypes$Cell_Type==ct), ])
  mean_nonpain_gex_nonct <- rowMeans(tpm_nonpain[, nontubes])
  test_result_ctVSothers <- t.test(mean_nonpain_gex_ct, mean_nonpain_gex_nonct, paired = T, alternative="greater")
  pval<-test_result_ctVSothers$p.value
  print(paste(ct, "vs Others Non-Pain Paired T test Alternative greater p-value:", pval))
  results <- append(results, -log(pval))
}
results[is.infinite(results)] <- NA
names(results) <- unique_cell_types

data_forplot <- data.frame(Category = names(results), Value = results)

p<- ggplot(data_forplot, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(title = "One vs Others Non-Pain Paired T test Alternative greater", x = "Cell Type", y = "-log p-value")

ggsave(plot = p, 
       filename=paste(save_dir, "One vs Others Non-Pain Paired T test Alternative greater_v3.png",sep=""),
       width = 12, height = 6, dpi = 180, units = "in", device='png')
