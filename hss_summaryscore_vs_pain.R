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



s <- top_k.Kendall
sampled_df <- df_ls_sorted[1:s, ]
genes <- row.names(sampled_df)

gex_zscores_sampled <- gex_zscores[genes,]

set.seed(42)
gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=TRUE) # # NOTE: perplexity needs to be consistent!

df_sampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                              y = pain_scores$pain_scores,
                              ps = pain_scores$pain_scores)

kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")

df_sampled.tsne <- data.frame(x = -sign(as.numeric(kendall.result$estimate[1]))*gex_zscores_sampled.tsne$Y,
                              y = pain_scores$pain_scores,
                              ps = pain_scores$pain_scores)

df_sampled.tsne.summary_score <- data.frame(gene_summary_score = -sign(as.numeric(kendall.result$estimate[1]))*gex_zscores_sampled.tsne$Y,
                                            pain_score = pain_scores$pain_scores,
                                            row.names = pain_scores$X)

#write.csv(df_sampled.tsne.summary_score,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/pain_related_top',s,'_gene_summary_scores.csv',sep=''))
# colnames(gex_zscores_sampled) <- pain_scores$X
# write.csv(gex_zscores_sampled,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/gex_zscore_pain_backup.csv',sep=''))


sp_ls.tsne<-ggplot(df_sampled.tsne, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE, colour="orange", size=7) + 
  #geom_point(aes(size=10), show.legend=FALSE, colour="black") + 
  geom_smooth(method = "lm", se = FALSE, colour="black") +
  theme(#text = element_text(size = 30),
        #axis.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size = 30),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))
#+
  #geom_text(hjust=0, vjust=-1.5) # The pain scores are for the y-axis. There is no need to repeat.

kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")

pa.kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")
pa.pearson.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="pearson")
pa.spearman.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="spearman")

figure_title <- paste('Kendall Correlation test P-value:',format(kendall.result$p.value,scientific=TRUE))
print("GbGMI result")
print(paste('Kendall Correlation test P-value:',kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',kendall.result$estimate))
gbgmi.p.value <- kendall.result$p.value
gbgmi.coef <- kendall.result$estimate

#+ggtitle(figure_title)
sp_ls.tsne + 
  #scale_x_continuous(name="815-Gene Subset Expression Summary Score") + 
  scale_y_continuous(limits=c(0, 60)) + #name="HOOS/KOOS Pain Score", 
  scale_x_continuous(limits=c(-150, 150)) +
  labs(x=NULL, y=NULL)

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_summary_scores_vs_pain_scores_pain_genes_Fig2",".png",sep=""), 
       width = 7, height = 10, dpi = 300, units = "in", device='png')
#dev.off()
######

set.seed(42)
gex_zscores.tsne <- Rtsne(t(gex_zscores), dims = 1, perplexity=perp, verbose=TRUE) # NOTE: perplexity needs to be consistent!

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                              y = pain_scores$pain_scores,
                              ps = pain_scores$pain_scores)

sp_ls.tsne<-ggplot(df.tsne, aes(x, y, label=ps)) +
  geom_point(aes(size=16), shape=18, show.legend=FALSE,colour="blue", size=8) +
  #geom_point(aes(size=10), show.legend=FALSE,colour="black") +
  geom_smooth(method = "lm", se = FALSE, color="black") +
  theme(#text = element_text(size = 30),
        #axis.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size = 30),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

all.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
all.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
all.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

print("All-gene result")
print(paste('Kendall Correlation test P-value:',all.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',all.kendall.result$estimate))

figure_title <- paste('Kendall Correlation test P-value:',format(kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne + 
  #scale_x_continuous(name="2,227-Gene Subset Expression Summary Score") + 
  scale_y_continuous(limits=c(0, 60)) +  #name="HOOS/KOOS Pain Score", 
  scale_x_continuous(limits=c(-150, 150)) +
  labs(x=NULL, y=NULL)
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_summary_scores_vs_pain_scores_low_infla_genes_paper_Fig2",".png",sep=""), 
       width = 7, height = 10, dpi = 300, units = "in", device='png')
#dev.off()

print("GbGMI result")
print(paste('Kendall Correlation test P-value:',gbgmi.p.value))
print(paste('Kendall Correlation coefficient:',gbgmi.coef))

print("Low inflammation genes result")
print(paste('Kendall Correlation test P-value:',all.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',all.kendall.result$estimate))




