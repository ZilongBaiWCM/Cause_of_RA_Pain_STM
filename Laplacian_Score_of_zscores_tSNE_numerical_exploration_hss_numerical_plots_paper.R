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

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_padj_<0.01log2FC_<0.csv',sep=""))
row.names(gex_zscores) <- gex_zscores$X
gex_zscores <- gex_zscores[,2:23]

pain_scores <- read.csv(paste(input_dir, 'pain_scores.csv', sep=""))

set.seed(42)
gex_zscores.umap = umap(t(gex_zscores))

df <- data.frame(x = gex_zscores.umap$layout[,1],
                 y = gex_zscores.umap$layout[,2],
                 ps = pain_scores$pain_scores)

sp2<-ggplot(df, aes(x, y, color=ps)) +
  geom_point(aes(size=ps))
sp2+scale_color_gradient(low="red", high="green")

# Load the filtered genes. Then do numerical exploration with parameter grid-search.
data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"
gex_ENSG_Symbols <- read.csv(paste(data_dir, 'gex_ENSG_Symbols.csv', sep=""))
row.names(gex_ENSG_Symbols) <- gex_ENSG_Symbols$X

ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
filter_condition <- '_zscores_bulk_padj_<0.01log2FC_<0'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores', filter_condition,'.csv', sep=""))
row.names(ls_result) <- ls_result$X

genes <- intersect(row.names(gex_ENSG_Symbols), row.names(ls_result))

df_ls <- merge(gex_ENSG_Symbols[genes,], ls_result[genes,])

attach(df_ls)
df_ls_sorted <- df_ls[order(Laplacian_Score),]
detach(df_ls)

row.names(df_ls_sorted) <- df_ls_sorted$X
df_ls_sorted <- df_ls_sorted[,2:4]

min_genes <- 1
max_genes <- dim(df_ls_sorted)[1]

df_embedded.vs.pains_cor <- data.frame()

perp <- 7

for(s in min_genes:max_genes){
  sampled_df <- df_ls_sorted[1:s, ]
  genes <- row.names(sampled_df) 
  gex_zscores_sampled <- gex_zscores[genes,]
  set.seed(42)
  gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=TRUE) #
  df_sampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                                y = pain_scores$pain_scores,
                                ps = pain_scores$pain_scores)
  pearson.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="pearson")
  kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")
  spearman.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="spearman")
  df_embedded.vs.pains_cor<-rbind(df_embedded.vs.pains_cor, data.frame(Pearson=pearson.result$p.value, Kendall=kendall.result$p.value, Spearman=spearman.result$p.value, combined = pearson.result$p.value + kendall.result$p.value )) #+ spearman.result$p.value
}

row.names(df_embedded.vs.pains_cor) <- min_genes:max_genes

top_k.Pearson <- row.names(df_embedded.vs.pains_cor)[which.min(df_embedded.vs.pains_cor$Pearson)]
top_k.Kendall <- row.names(df_embedded.vs.pains_cor)[which.min(df_embedded.vs.pains_cor$Kendall)]
top_k.Spearman <- row.names(df_embedded.vs.pains_cor)[which.min(df_embedded.vs.pains_cor$Spearman)]

top_k.combined <- row.names(df_embedded.vs.pains_cor)[which.min(df_embedded.vs.pains_cor$combined)]

# Plot the top k genes 
xind <- 1:dim(df_embedded.vs.pains_cor)[1]
ggplot() + 
  geom_segment(aes(x = xind, xend = xind, y = 0, yend = -log(Kendall)), data = df_embedded.vs.pains_cor) +
  geom_point(aes( x = xind, y = -log(Kendall)), colour = "orange", size = 3.5, data = df_embedded.vs.pains_cor) + 
  #geom_point(aes( x = xind, y = -log(Kendall)), colour = "black", size = 3.5, data = df_embedded.vs.pains_cor) + 
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))+
  #xlab("k (top k-gene subset)") +
  #ylab("-log(P)")+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(0,7))

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Kendall_test_p_value_vs_top_k_summary_score_LS_sorted_HSS_Fig2D",".png",sep=""), width = 6, height = 6, dpi = 300, units = "in", device='png')

df_single_gene.vs.pains_cor <- data.frame()

for(s in min_genes:max_genes){
  #sampled_df <- rownames(df_ls_sorted)[s]
  #gene <- sampled_df
  gene <- rownames(df_ls_sorted)[s]
  gex_zscores_sampled <- gex_zscores[gene,]
  gex_zscores_sampled <- t(gex_zscores_sampled[1,])
  df_single_gene <- data.frame(x = as.vector(t(gex_zscores_sampled)),
                                y = pain_scores$pain_scores,
                                ps = pain_scores$pain_scores)
  pearson.result <- cor.test(df_single_gene$x, df_single_gene$y, method="pearson")
  kendall.result <- cor.test(df_single_gene$x, df_single_gene$y, method="kendall")
  spearman.result <- cor.test(df_single_gene$x, df_single_gene$y, method="spearman")
  df_single_gene.vs.pains_cor<-rbind(df_single_gene.vs.pains_cor, data.frame(Pearson=pearson.result$p.value, Kendall=kendall.result$p.value, Spearman=spearman.result$p.value, combined = pearson.result$p.value + kendall.result$p.value )) #+ spearman.result$p.value
}

# Plot the k-th gene 
xind <- 1:dim(df_single_gene.vs.pains_cor)[1]
ggplot() + 
  geom_segment(aes(x = xind, xend = xind, y = 0, yend = -log(Kendall)), data = df_single_gene.vs.pains_cor) +
  geom_point(aes( x = xind, y = -log(Kendall)), colour = "pink", size = 2, data = df_single_gene.vs.pains_cor) + 
  theme(text = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.y = element_text(size = 30))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("k-th gene (sorted by Laplacian Scores in ascending order)") +
  ylab("log(P-value in Kendall Correlation test between single gene expression and pain scores)") +
  scale_y_continuous(limits=c(0,7))

#ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Kendall_test_p_value_vs_kth_single_gene_LS_sorted_HSS_Fig2BK",".png",sep=""), width = 6, height = 6, dpi = 300, units = "in", device='png')


# Plot the k-th gene. Laplacian Scores
xind <- 1:dim(df_ls_sorted)[1]
ggplot() + 
  #geom_segment(aes(x = xind, xend = xind, y = 0, yend = Laplacian_Score), data = df_ls_sorted) +
  #geom_point(aes( x = xind, y = Laplacian_Score, colour = Laplacian_Score), size = 3, data = df_ls_sorted,show.legend = FALSE) + 
  geom_point(aes( x = xind, y = Laplacian_Score, colour = "darkblue"), size = 3, data = df_ls_sorted,show.legend = FALSE, color="darkblue") + 
  #geom_point(aes( x = xind, y = Laplacian_Score, colour = "darkblue"), size = 3, data = df_ls_sorted,show.legend = FALSE, color="black") + 
  # scale_colour_gradient2(
  #   low = "gold",
  #   mid = "white",
  #   high = "blue",
  #   midpoint = (max(df_ls_sorted$Laplacian_Score)+min(df_ls_sorted$Laplacian_Score))/2)+
  theme(#text = element_text(size=rel(0.5)),
    #axis.text = element_text(size=120),
    #axis.text.x = element_text(size=20),
    axis.title.x = element_blank(),
    #axis.text.y = element_text(size=20),
    axis.title.y = element_blank()
  ) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))+
  # xlab("k (gene index sorted by Laplacian Scores in ascending order)") +
  # ylab("Laplacian Score") +
  labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(0.7,1))

ggsave(file=paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Laplacian_Score_kth_single_gene_LS_sorted_HSS_Fig2C",".png",sep=""), width = 6, height = 6, dpi = 300, units = "in", device='png')


