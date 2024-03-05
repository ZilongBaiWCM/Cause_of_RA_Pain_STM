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

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_vs_nuc_all_patients_5000_genes.csv',sep=""))
row.names(gex_zscores) <- gex_zscores$X
gex_zscores <- gex_zscores[,2:dim(gex_zscores)[2]]

nuc_density_scores <- read.csv(paste(input_dir, 'nuc_density_scores_all_patients_5000_genes.csv', sep=""))

set.seed(42)
gex_zscores.umap = umap(t(gex_zscores))

df <- data.frame(x = gex_zscores.umap$layout[,1],
                 y = gex_zscores.umap$layout[,2],
                 nds = nuc_density_scores$nuc_density_scores)

sp2<-ggplot(df, aes(x, y, color=nds)) +
  geom_point(aes(size=nds))
sp2+scale_color_gradient(low="red", high="green")

# Load the filtered genes. Then do numerical exploration with parameter grid-search.
data_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"
gex_ENSG_Symbols <- read.csv(paste(data_dir, 'gex_ENSG_Symbols.csv', sep=""))
row.names(gex_ENSG_Symbols) <- gex_ENSG_Symbols$X

ls_dir <- '/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/Graph_Regularized_Feature_Selection/'
ls_result<- read.csv(paste(ls_dir,'Laplacian_Scores_bulk_vs_nuc_all_patients_5000_genes.csv', sep=""))
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

df_embedded.vs.nuc_cor <- data.frame()

perp <- 2 # mean degree of the weighted graph for nuc is: 2.067431536475631

for(s in min_genes:max_genes){
  sampled_df <- df_ls_sorted[1:s, ]
  genes <- row.names(sampled_df) 
  gex_zscores_sampled <- gex_zscores[genes,]
  set.seed(42)
  gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=TRUE) #
  df_sampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                                y = nuc_density_scores$nuc_density_scores,
                                nds = nuc_density_scores$nuc_density_scores)
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

#NOTE: for tsne, perplexity=6, dims=1 projection 
# Gives top_k.Kendall = 585. This gene list has neuron projection enriched. 

s <- top_k.Kendall
sampled_df <- df_ls_sorted[1:s, ]
genes <- row.names(sampled_df)

gex_zscores_sampled <- gex_zscores[genes,]

set.seed(42)
gex_zscores_sampled.tsne <- Rtsne(t(gex_zscores_sampled), dims = 1, perplexity=perp, verbose=TRUE) # # NOTE: perplexity needs to be consistent!

df_sampled.tsne <- data.frame(x = gex_zscores_sampled.tsne$Y,
                              y = nuc_density_scores$nuc_density_scores,
                              nds = nuc_density_scores$nuc_density_scores)

kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")

if (sign(as.numeric(kendall.result$estimate[1])) > 0){
  sign_sampled <- 1
}else{
  sign_sampled <- -1
}

df_sampled.tsne <- data.frame(x = sign_sampled*gex_zscores_sampled.tsne$Y,
                              y = nuc_density_scores$nuc_density_scores,
                              nds = nuc_density_scores$nuc_density_scores)

df_sampled.tsne.summary_score <- data.frame(gene_summary_score = sign_sampled*gex_zscores_sampled.tsne$Y,
                                            nuc_density_scores = nuc_density_scores$nuc_density_scores,
                                            row.names = nuc_density_scores$X)

#write.csv(df_sampled.tsne.summary_score,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/pain_related_top',s,'_gene_summary_scores.csv',sep=''))

sp_ls.tsne<-ggplot(df_sampled.tsne, aes(x, y, label=nds)) +
  geom_point(aes(size=10), show.legend=FALSE, colour="orange", size = 7) + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))

#+
#geom_text(hjust=0, vjust=-1.5) # The pain scores are for the y-axis. There is no need to repeat.

kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")

nds.kendall.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="kendall")
nds.pearson.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="pearson")
nds.spearman.result <- cor.test(df_sampled.tsne$x, df_sampled.tsne$y, method="spearman")

figure_title <- paste('Kendall Correlation test P-value:',format(nds.kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne+
  labs(x=NULL, y=NULL)+ # + scale_x_continuous(name="2713-Gene Subset Expression Summary Score") + 
  scale_y_continuous(limits=c(1500, 7000)) #name="Cell Density", 
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Fig2G_summary_scores_vs_average_cell_density",".png",sep=""), width = 5, height = 6, dpi = 300, units = "in", device='png')
#dev.off()
######

set.seed(42)
gex_zscores.tsne <- Rtsne(t(gex_zscores), dims = 1, perplexity=perp, verbose=TRUE) # NOTE: perplexity needs to be consistent!

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                      y = nuc_density_scores$nuc_density_scores,
                      nds = nuc_density_scores$nuc_density_scores)

kendall.result.all <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

if (sign(as.numeric(kendall.result.all$estimate[1])) > 0){
  sign_all <- 1
}else{
  sign_all <- -1
}

df.tsne <- data.frame(x = sign_all*gex_zscores.tsne$Y,
                              y = nuc_density_scores$nuc_density_scores,
                              nds = nuc_density_scores$nuc_density_scores)

sp_ls.tsne<-ggplot(df.tsne, aes(x, y, label=nds)) +
  geom_point(aes(size=10), show.legend=FALSE,colour="blue", size = 7) + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

pa.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
pa.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
pa.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

figure_title <- paste('Kendall Correlation test P-value:',format(kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne + 
  labs(x=NULL, y=NULL)+
  #scale_x_continuous(name="5000-Gene Subset Expression Summary Score") + 
  scale_y_continuous(limits=c(1500, 7000)) #name="Cell Density", 
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Fig2E_summary_scores_vs_average_cell_density_5000_genes",".png",sep=""), width = 5, height = 6, dpi = 300, units = "in", device='png')
#dev.off()



