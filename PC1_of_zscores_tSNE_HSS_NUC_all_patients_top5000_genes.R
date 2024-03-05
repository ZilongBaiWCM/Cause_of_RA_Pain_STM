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

# NOTE: This code is ongoing. In order to demonstrate the power of our method, we should not have filtered the genes to be the ones upregulated in DEA of high vs relatively low.

input_dir <- "/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/InputData/"

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_vs_nuc_all_patients_5000_genes.csv',sep=""))
row.names(gex_zscores) <- gex_zscores$X
gex_zscores <- gex_zscores[,2:dim(gex_zscores)[2]]

gex_zscores.pca <- prcomp(gex_zscores)

gex_zscores.pca.pc1 <- gex_zscores.pca$rotation[,'PC1']

nuc_density_scores <- read.csv(paste(input_dir, 'nuc_density_scores_all_patients_5000_genes.csv', sep=""))

df_sampled.pc1 <- data.frame(x = gex_zscores.pca.pc1,
                              y = nuc_density_scores$nuc_density_scores,
                              nds = nuc_density_scores$nuc_density_scores)

kendall.result <- cor.test(df_sampled.pc1$x, df_sampled.pc1$y, method="kendall")

df_sampled.pc1 <- data.frame(x = sign(as.numeric(kendall.result$estimate[1]))*gex_zscores.pca.pc1, #sign(as.numeric(kendall.result$estimate[1]))*
                             y = nuc_density_scores$nuc_density_scores,
                             nds = nuc_density_scores$nuc_density_scores)

df_sampled.pc1.summary_score <- data.frame(gex_zscores.pca.pc1 = sign(as.numeric(kendall.result$estimate[1]))*gex_zscores.pca.pc1, #sign(as.numeric(kendall.result$estimate[1]))*
                                            nuc_density_scores = nuc_density_scores$nuc_density_scores,
                                            row.names = nuc_density_scores$X)

#write.csv(df_sampled.pc1.summary_score,paste('/Users/baizilong/Documents/Dana-dataset/Exploration-results/Paper_Gene-Expression_Low-Inflammatory-Synovium_Pain/gex_pc1.vs.nuc_benchmarking_all_patients_5000_genes.csv',sep=''))

kendall.result <- cor.test(df_sampled.pc1$x, df_sampled.pc1$y, method="kendall")
pearson.result <- cor.test(df_sampled.pc1$x, df_sampled.pc1$y, method="pearson")
spearman.result <- cor.test(df_sampled.pc1$x, df_sampled.pc1$y, method="spearman")

#figure_title <- paste('Kendall Correlation test P-value:',format(kendall.result$p.value,scientific=TRUE))
figure_title <- paste('Kendall:',format(kendall.result$p.value, scientific=TRUE), 'Pearson:', format(pearson.result$p.value, scientific=TRUE)) #,scientific=TRUE

sp_pc1.pca<-ggplot(df_sampled.pc1, aes(x, y, color=nds, label=nds)) +
  geom_point(aes(size=10), show.legend=FALSE,colour="purple", size = 7) + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme(#text = element_text(size=rel(0.5)),
        #axis.text = element_text(size=120),
        #axis.text.x = element_text(size=20),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size=20),
        axis.title.y = element_blank())+
  #ggtitle(figure_title)+
  #scale_color_gradient(low="red", high="blue") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))+
  #scale_x_continuous(name="5000-Gene Subset Expression PC1") + 
  labs(x=NULL, y=NULL)+
  scale_y_continuous(limits=c(1500, 7000)) #name="Cell Density", 
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/Fig2F_pc1_vs_average_cell_density_all_patients_5000_genes",".png",sep=""), width = 5, height = 6, dpi = 300, units = "in", device='png')
