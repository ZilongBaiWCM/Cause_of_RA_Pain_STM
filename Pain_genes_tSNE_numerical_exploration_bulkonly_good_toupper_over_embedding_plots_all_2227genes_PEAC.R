# This script aims to find pathways that are "similar" to the Neuroactive Ligand-Receptor Interaction.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
rm(list=ls())

#library(SingleCellExperiment)
library(ggplot2)
library(magrittr)
library(DESeq2)
library(dplyr)
#install.packages("factoextra")
library(factoextra)
#BiocManager::install("clusterProfiler")
library(enrichplot)
library(org.Hs.eg.db) # Gene Ontology
library(clusterProfiler)
#BiocManager::install('M3C')
library(M3C) # TSNE
library(umap) # umap
library(cluster)
#install.packages('amap')
library(amap) # Kmeans
library(tidyverse)  # data manipulation
library(dendextend) # for comparing two dendrograms
# load library
library(ggrepel) # Volcano plot
library(tidyverse)
library(readxl)
#install.packages('xlsx')
library(xlsx)
library("readxl")

library(ggplot2)

library(Rtsne)

# Loading 815 pain-associated genes identified by my algorithm.
ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores', filter_condition,'.csv', sep=""))
row.names(ls_result) <- ls_result$X

# Loading PEAC gene expression data.
peac_dir <- '/Users/baizilong/Documents/Dana-dataset/Validation/'
exp_fn <- 'batch_corrected_unfiltered_with_gene_log_cpm_peac_04192022.txt'
exp_data <- read.table(paste(peac_dir, exp_fn, sep=""), header=T, row.names =1)

row.names(exp_data) <- gsub("\\..*","",exp_data$gene)

# Extract PEAC gene expression for pain-associated genes.
pain_exp_data <- exp_data[toupper(row.names(exp_data)) %in% toupper(row.names(ls_result)),]
print(dim(pain_exp_data))
if (length(unique(pain_exp_data$SYMBOL)) == dim(pain_exp_data)[1]){
  row.names(pain_exp_data) <- pain_exp_data$SYMBOL
}

# Loading PEAC meta data
meta_fn <- 'PEAC FASTQ METADATA.xlsx'
meta_data <- read_excel(paste(peac_dir, meta_fn, sep=""), sheet = 2)
# Focus on "fibroid"

vas_scores <- meta_data[c('Technology','Characteristics[VAS]', 'Characteristics[pathotype]')]

colnames(vas_scores) <- c('Technology', 'VAS', 'Pathotype')

# Merge fibroid and ungraded into low inflammatory subgroup
vas_scores['Inflammation.Level'] <- 'NaN'
for (r in row.names(vas_scores)){
  if (vas_scores[r, 'Pathotype'] == 'lymphoid'){
    vas_scores[r, 'Inflammation.Level'] <- 'High: lymphoid'
  }
  if (vas_scores[r, 'Pathotype'] == 'myeloid'){
    vas_scores[r, 'Inflammation.Level'] <- 'Mid: myeloid'
  }
  if (vas_scores[r, 'Pathotype'] == 'fibroid' | vas_scores[r, 'Pathotype'] == 'ungraded' ){
    vas_scores[r, 'Inflammation.Level'] <- 'Low: fibroid+ungraded'
  }
}

pathotype <- 'all'

vas_scores <- unique(vas_scores) # Only keep unique rows.
vas_scores <- as.data.frame(vas_scores)
names.vas_scores <- vas_scores$Technology
vas_scores_only <- vas_scores[c('VAS')]
row.names(vas_scores_only) <- names.vas_scores

write.csv(vas_scores_only, paste(peac_dir, pathotype, '_vas_scores.csv', sep=""))

# Extract PEAC gene expression for our pain-associated genes and fibroid samples
samples <- intersect(colnames(pain_exp_data), names.vas_scores)
pain_exp_data_sampled <- pain_exp_data[, c('SYMBOL.1', samples)]
row.names(pain_exp_data_sampled) <- pain_exp_data_sampled$SYMBOL
pain_gex_sampled <- subset(pain_exp_data_sampled, select=-c(SYMBOL.1))

# Z-scoring gene expression
gex_zscores <- apply(pain_gex_sampled, MARGIN = 1, FUN = scale )
gex_zscores <- t(gex_zscores)

gex_zscores_sampled <- gex_zscores
colnames(gex_zscores_sampled) <- colnames(pain_gex_sampled)

rownames(vas_scores) <- rownames(vas_scores)
colnames(gex_zscores_sampled) <- rownames(vas_scores)

perp = 20 # Mean degree: 20.04210043386775  # all

set.seed(42)
gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=TRUE) # # NOTE: perplexity needs to be consistent!

df_sampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                              y = vas_scores$VAS,
                              ps = vas_scores$VAS,
                              inflammation = vas_scores$Inflammation.Level)

unique_inflammation <- unique(df_sampled.tsne$inflammation)

coloring <- c('purple', 'green', 'blue')
names(coloring) <- unique_inflammation

for(infla in unique_inflammation){
  df_subsampled.tsne <- df_sampled.tsne[df_sampled.tsne$inflammation == infla,]
  
  kendall.result <- cor.test(df_subsampled.tsne$x, df_subsampled.tsne$y, method="kendall")
  pearson.result <- cor.test(df_subsampled.tsne$x, df_subsampled.tsne$y, method="pearson")
  spearman.result <- cor.test(df_subsampled.tsne$x, df_subsampled.tsne$y, method="spearman")
  
  print(infla)
  kendall_output <- paste('kendall:', 'tau=', kendall.result$estimate, 'p-value=', kendall.result$p.value)
  
  pearson_output <- paste('pearson:', 'cor=', pearson.result$estimate, 'p-value=', pearson.result$p.value)
  
  spearman_output <- paste('spearman:', 'rho=', kendall.result$estimate, 'p-value=', kendall.result$p.value)
  
  print(kendall_output)
  print(pearson_output)
  print(spearman_output)
  
  sign_info <- pearson.result$estimate[1]
  if (sign_info < 0 ){
    sign_info <- -1
  }else{
    sign_info <- 1
  }
  
  # sign_info <- pearson.result$estimate[1]
  # if(sign_info == 0){
  #   print('sign_info==0. Not changing sign. Correlation estimate == 0')
  #   sign_info <- 1
  # }else{
  #   sign_info <- sign(as.numeric(sign_info))
  # }
  
  df_subsampled.tsne.signed <- data.frame(x = sign_info*df_subsampled.tsne$x,
                                y = df_subsampled.tsne$y,
                                ps = df_subsampled.tsne$ps)
  
  kendall.result.signed <- cor.test(df_subsampled.tsne.signed$x, df_subsampled.tsne.signed$y, method="kendall")
  pearson.result.signed <- cor.test(df_subsampled.tsne.signed$x, df_subsampled.tsne.signed$y, method="pearson")
  spearman.result.signed <- cor.test(df_subsampled.tsne.signed$x, df_subsampled.tsne.signed$y, method="spearman")
  
  print(infla)
  kendall_output.signed <- paste('Kendall:', 'tau=', format(kendall.result.signed$estimate,scientific=F), 'p-value=', format(kendall.result.signed$p.value,scientific=F))
  
  pearson_output.signed <- paste('Pearson:', 'cor=', format(pearson.result.signed$estimate,scientific=F), 'p-value=', format(pearson.result.signed$p.value,scientific=F))
  
  spearman_output.signed <- paste('Spearman:', 'rho=', format(kendall.result.signed$estimate,scientific=F), 'p-value=', format(kendall.result.signed$p.value,scientific=F))
  
  df_subsampled.tsne.summary_score.signed <- data.frame(gene_summary_score = sign_info*df_subsampled.tsne$x,
                                              pain_score = df_subsampled.tsne$y,
                                              row.names = row.names(df_subsampled.tsne))
  
  print(kendall_output.signed)
  print(pearson_output.signed)
  print(spearman_output.signed)
  
  figure_title <- paste('Summary Score vs VAS in', infla)
  x_label = paste('Relatively Low-inflammatory ',dim(pain_exp_data)[1],'-Gene Expression Summary Score')
  sp_ls.tsne<-ggplot(df_subsampled.tsne.signed, aes(x, y)) +
    geom_point(aes(size=15), shape=18, show.legend=FALSE, color=coloring[infla], size=7) + 
    geom_smooth(method = "lm", se = FALSE) +
    theme(#text = element_text(size=rel(0.5)),
      #axis.text = element_text(size=120),
      #axis.text.x = element_text(size=20),
      axis.title.x = element_blank(),
      #axis.text.y = element_text(size=20),
      axis.title.y = element_blank()
    ) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       axis.text = element_text(size=30))+
    #ggtitle(figure_title)+
    #scale_color_manual(values = coloring)
    #scale_x_continuous(name=x_label) + 
    scale_y_continuous(limits=c(0, 100))+
    labs(x=NULL, y=NULL)
  
  #ggsave(plot = sp_ls.tsne, paste("/Users/baizilong/Documents/Dana-dataset/RA_paper/peac_validation_lowinf_pain_genes_summary_scores_vs_VAS_",infla,"_Fig2.png",sep=""))
  
}

