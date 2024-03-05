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
s = 815
pain_genes <- read.csv(paste(ls_dir,'Laplacian_Scores_of_zscores_bulkonly_top',s,'.csv',sep=""), row.names = 1)
print(dim(pain_genes))

# Set gex and pain score directory.
data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"

unique_inflammation <- c('highinf', 'lowinf','all') #c('High', 'Low (combined)','all')
#
inflatypes <- c('highinf', 'lowinf','all')
names(inflatypes) <- unique_inflammation
#
coloring <- c('#5A5A5A', 'orange', 'blue')
names(coloring) <- unique_inflammation
#
#ud = 3.1 #3.1, 6.1
degrees <- c(3,7,11) 
names(degrees) <- unique_inflammation


for(infla in unique_inflammation){
  inflat = inflatypes[infla]
  perp = degrees[infla]
  
  print(infla)
  print(inflat)
  print(perp)
  # Load z-scored gene expression.
  gex_zscores <- read.csv(paste(data_dir, 'GEX_zscores_bulk_pain815_', inflat, '_HSS_RA_VAS.csv', sep=""), row.names = 1)
  
  # Load pain scores
  vas_scores <- read.csv(paste(data_dir, 'pain_scores_hss_ra_vas_', inflat, '.csv', sep=""), row.names = 1)
  
  samples <- intersect(colnames(gex_zscores), rownames(vas_scores))
  genes <- intersect(rownames(gex_zscores), rownames(pain_genes))
  
  gex_zscores_sampled <- gex_zscores[genes, samples]
  
  print(rownames(vas_scores) == colnames(gex_zscores_sampled))
  
  set.seed(42)
  gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=F) # # NOTE: perplexity needs to be consistent!
  
  df_subsampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                                y = vas_scores$pain_scores,
                                ps = vas_scores$pain_scores,
                                inflammation = vas_scores$inflammation_level)
  
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
  
  figure_title <- paste('Summary Score vs VAS (HSS) in', infla)
  sp_ls.tsne<-ggplot(df_subsampled.tsne.signed, aes(x, y)) +
    geom_point(aes(size=15), show.legend=FALSE, color=coloring[infla], size=7) + 
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
    #scale_x_continuous(name="Pain-associated 738-Gene Expression Summary Score") + 
    scale_y_continuous(limits=c(0, 10.1))+
    labs(x=NULL, y=NULL)
  
  if (infla=='lowinf'){
    sp_ls.tsne <- sp_ls.tsne + scale_x_continuous(limits=c(-150, 150))
  }
  
  if (infla=='all'){
    sp_ls.tsne <- sp_ls.tsne + scale_x_continuous(limits=c(-800, 800))
  }
  
  print(range(df_subsampled.tsne.signed$x))
  
  ggsave(plot = sp_ls.tsne, paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_validation_pain_genes_summary_scores_vs_HSS_RA_VAS_",infla,"_Fig2_extended_v2.png",sep=""),
         width = 7, height = 10, dpi = 150, units = "in", device='png')
  
}

