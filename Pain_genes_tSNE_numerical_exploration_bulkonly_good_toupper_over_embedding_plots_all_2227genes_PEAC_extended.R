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

# Loading 2227 low inflammatory genes identified by my algorithm.
ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores', filter_condition,'.csv', sep=""))
row.names(ls_result) <- ls_result$X



#set.seed(42)
# Set gex and pain score directory.
data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"

unique_inflammation <- c('High: Lymphoid', 'Mixed: Myeloid','Low: Fibroid+Ungraded','All', 'Lymphoid+Fibroid+Ungraded', 'Low:Myeloid+Fibroid', 'Low:Myeloid+Fibroid+Ungraded')
#
pathotypes <- c('lymphoid', 'myeloid','fibroidungraded','all', 'lymfibun', 'myefib', 'myefibun')
names(pathotypes) <- unique_inflammation
#
coloring <- c('blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue')
names(coloring) <- unique_inflammation
#
# Pathotype: lymphoid
# Graph mean degree: 11.974016523020804
# Pathotype: myeloid
# Graph mean degree: 4.637079703790972
# Pathotype: fibroidungraded
# Graph mean degree: 4.869236350279577
# Pathotype: all
# Graph mean degree: 20.04210043386775
# Pathotype: lymfibun
# Graph mean degree: 16.385663571722468
# Pathotype: myefib
# Graph mean degree: 7.629047346168079

#ud = 3
degrees <- c(11,5,3,17,15,8,7) 
names(degrees) <- unique_inflammation

for(infla in unique_inflammation){
  # 
  patho = pathotypes[infla]
  perp = degrees[infla]
  
  print(infla)
  print(patho)
  print(perp)
  # Load z-scored gene expression.
  gex_zscores <- read.csv(paste(data_dir, 'GEX_zscores_bulk_', patho, '_PEAC_v2.csv', sep=""), row.names = 1)
  
  # Load pain scores
  vas_scores <- read.csv(paste(data_dir, 'pain_scores_peac_', patho, '.csv', sep=""), row.names = 1)
  
  samples <- intersect(colnames(gex_zscores), rownames(vas_scores))
  genes <- intersect(rownames(gex_zscores), rownames(ls_result))
  
  gex_zscores_sampled <- gex_zscores[genes, samples]
  vas_scores <- vas_scores[samples,]
  
  #print(paste('gex dimentionality:', dim(gex_zscores_sampled)))
  # print(colnames(gex_zscores_sampled))
  # print(rownames(gex_zscores_sampled))
  
  set.seed(42)

  gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=F) # # NOTE: perplexity needs to be consistent!
  
  df_subsampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                                   y = vas_scores$VAS,
                                   ps = vas_scores$VAS,
                                   inflammation = vas_scores$Pathotype)
  
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
  
  # df_subsampled.tsne.summary_score.signed <- data.frame(gene_summary_score = sign_info*df_subsampled.tsne$x,
  #                                             pain_score = df_subsampled.tsne$y,
  #                                             row.names = row.names(df_subsampled.tsne))
  
  print(kendall_output.signed)
  print(pearson_output.signed)
  print(spearman_output.signed)
  
  print(range(df_subsampled.tsne.signed$x))
  
  figure_title <- paste('Summary Score vs VAS in', infla)
  #x_label = paste('Relatively Low-inflammatory ',dim(pain_exp_data)[1],'-Gene Expression Summary Score')
  sp_ls.tsne<-ggplot(df_subsampled.tsne.signed, aes(x, y)) +
    geom_point(aes(size=16), shape=18, show.legend=FALSE, color=coloring[infla], size=8) + 
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
    scale_x_continuous(limits=c(-610, 610))+
    labs(x=NULL, y=NULL)
  
  ggsave(plot = sp_ls.tsne, paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/peac_validation_lowinf_pain_genes_summary_scores_vs_VAS_",infla,"_Fig2_extended.png",sep=""),
         width = 7, height = 10, dpi = 150, units = "in", device='png')
  
}

