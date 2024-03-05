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

library(ggplot2)

library(Rtsne)

input_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_vs_nuc_all_patients_5000_genes.csv',sep=""), row.names = 1)

nuc_density <- read.csv(paste(input_dir, 'nuc_density_scores_all_patients_5000_genes.csv', sep=""))

set.seed(42)
gex_zscores.umap = umap(t(gex_zscores))

df <- data.frame(x = gex_zscores.umap$layout[,1],
                 y = gex_zscores.umap$layout[,2],
                 ps = nuc_density$nuc_density_scores)

sp2<-ggplot(df, aes(x, y, color=ps)) +
  geom_point(aes(size=ps))
sp2+scale_color_gradient(low="red", high="green")

# Load the filtered genes. Then do numerical exploration with parameter grid-search.
data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"
gex_ENSG_Symbols <- read.csv(paste(data_dir, 'gex_ENSG_Symbols.csv', sep=""),row.names=1)
#row.names(gex_ENSG_Symbols) <- gex_ENSG_Symbols$X

ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores_bulk_vs_nuc_all_patients_5000_genes.csv', sep=""),row.names=1)
#row.names(ls_result) <- ls_result$X

genes <- intersect(row.names(gex_ENSG_Symbols), row.names(ls_result))

df_ls <- merge(gex_ENSG_Symbols[genes,], ls_result[genes,,drop=F],by = 'row.names')

attach(df_ls)
df_ls_sorted <- df_ls[order(Laplacian_Score),]
detach(df_ls)

row.names(df_ls_sorted) <- df_ls_sorted$Row.names
df_ls_sorted <- df_ls_sorted[,2:4]

min_genes <- 1
max_genes <- dim(df_ls_sorted)[1]

df_embedded.vs.nuc_cor <- data.frame()

perp <- 2 # average degree:2.067431536475631

for(s in min_genes:max_genes){
  sampled_df <- df_ls_sorted[1:s, ]
  genes <- row.names(sampled_df) 
  gex_zscores_sampled <- gex_zscores[genes,]
  set.seed(42)
  gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=TRUE) #
  df_sampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                                y = nuc_density$nuc_density_scores,
                                ps = nuc_density$nuc_density_scores)
  pearson.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="pearson")
  kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")
  spearman.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="spearman")
  df_embedded.vs.nuc_cor<-rbind(df_embedded.vs.nuc_cor, data.frame(Pearson=pearson.result$p.value, Kendall=kendall.result$p.value, Spearman=spearman.result$p.value, combined = pearson.result$p.value + kendall.result$p.value )) #+ spearman.result$p.value
}

row.names(df_embedded.vs.nuc_cor) <- min_genes:max_genes

top_k.Pearson <- row.names(df_embedded.vs.nuc_cor)[which.min(df_embedded.vs.nuc_cor$Pearson)]
top_k.Kendall <- row.names(df_embedded.vs.nuc_cor)[which.min(df_embedded.vs.nuc_cor$Kendall)]
top_k.Spearman <- row.names(df_embedded.vs.nuc_cor)[which.min(df_embedded.vs.nuc_cor$Spearman)]

top_k.combined <- row.names(df_embedded.vs.nuc_cor)[which.min(df_embedded.vs.nuc_cor$combined)]

# Plot the top k genes 
xind <- 1:dim(df_embedded.vs.nuc_cor)[1]
ggplot() + 
  geom_segment(aes(x = xind, xend = xind, y = 0, yend = -log(Kendall)), data = df_embedded.vs.nuc_cor) +
  geom_point(aes( x = xind, y = -log(Kendall)), colour = "orange", size = 3.5, data = df_embedded.vs.nuc_cor) + 
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))+
  # xlab("Top k genes (in ascending Laplacian Score)") +
  # ylab("-log(P-value)")+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(0,10))

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Fig2D_Kendall_test_p_value_vs_top_k_summary_score_LS_sorted_NUC",".png",sep=""), width = 8, height = 6, dpi = 300, units = "in", device='png')

df_single_gene.vs.nuc_cor <- data.frame()

for(s in min_genes:max_genes){
  #sampled_df <- rownames(df_ls_sorted)[s]
  #gene <- sampled_df
  gene <- rownames(df_ls_sorted)[s]
  gex_zscores_sampled <- gex_zscores[gene,]
  gex_zscores_sampled <- t(gex_zscores_sampled[1,])
  df_single_gene <- data.frame(x = as.vector(t(gex_zscores_sampled)),
                                y = nuc_density$nuc_density_scores,
                                ps = nuc_density$nuc_density_scores)
  pearson.result <- cor.test(df_single_gene$x, df_single_gene$y, method="pearson")
  kendall.result <- cor.test(df_single_gene$x, df_single_gene$y, method="kendall")
  spearman.result <- cor.test(df_single_gene$x, df_single_gene$y, method="spearman")
  df_single_gene.vs.nuc_cor<-rbind(df_single_gene.vs.nuc_cor, data.frame(Pearson=pearson.result$p.value, Kendall=kendall.result$p.value, Spearman=spearman.result$p.value, combined = pearson.result$p.value + kendall.result$p.value )) #+ spearman.result$p.value
}

# Plot the k-th gene 
xind <- 1:dim(df_single_gene.vs.nuc_cor)[1]
ggplot() + 
  geom_segment(aes(x = xind, xend = xind, y = 0, yend = -log(Kendall)), data = df_single_gene.vs.nuc_cor) +
  geom_point(aes( x = xind, y = -log(Kendall)), colour = "pink", size = 3.5, data = df_single_gene.vs.nuc_cor) +
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))+
  # xlab("k-th gene (in ascending Laplacian Score)") +
  # ylab("-log(P-value)") +
  labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(0,10))

#ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Kendall_test_p_value_vs_kth_single_gene_LS_sorted_NUCBK",".png",sep=""), width = 8, height = 6, dpi = 300, units = "in", device='png')


# Plot the k-th gene. Laplacian Scores
xind <- 1:dim(df_ls_sorted)[1]
ggplot() + 
  #geom_segment(aes(x = xind, xend = xind, y = 0, yend = Laplacian_Score), data = df_ls_sorted) +
  geom_point(aes( x = xind, y = Laplacian_Score), colour = "blue", size = 3.5, shape=10, data = df_ls_sorted) + 
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))+
  # xlab("k-th gene (in ascending Laplacian Score)") +
  # ylab("Laplacian Score") +
  labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,5200))

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Fig2C_Laplacian_Score_kth_single_gene_LS_sorted_NUC",".png",sep=""), width = 6, height = 6, dpi = 300, units = "in", device='png')


