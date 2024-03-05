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

pain_scores <- read.csv(paste(input_dir, 'pain_scores_hooskoos.csv', sep=""), row.names = 1)

# All patients

gex <- read.csv(paste(input_dir, 'GEX_TF_bulk.csv',sep=""), row.names = 1)

gex <- data.frame(t(gex))

pain_scores_all <- pain_scores[rownames(gex),]

df <- data.frame(x = gex$TF,
                y = pain_scores_all$pain_scores,
                ps = pain_scores_all$pain_scores)

all.kendall.result <- cor.test(df$x, df$y, method="kendall")
all.pearson.result <- cor.test(df$x, df$y, method="pearson")
all.spearman.result <- cor.test(df$x, df$y, method="spearman")

tf_vs_pain<-ggplot(df, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="blue", size=7) +
  #geom_point(aes(size=10), show.legend=FALSE,colour="black") +
  geom_smooth(method = "lm", se = FALSE, color="black") +
  theme(#text = element_text(size = 30),
        #axis.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        #axis.text.y = element_text(size = 30),
        axis.title.y = element_blank())+
  labs(x=NULL, y=NULL)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))



print("All-patient HOOS/KOOS result")

print(paste('Kendall Correlation test P-value:',all.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',all.kendall.result$estimate))

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_tf_expression_vs_pain_scores_allinf",".png",sep=""),
       width = 7, height = 10, dpi = 300, units = "in", device='png')
#dev.off()

#==============================================================================================================================#

# Low inflammatory patients

gex <- read.csv(paste(input_dir, 'GEX_TF_bulk_lowinf.csv',sep=""), row.names = 1)

gex <- data.frame(t(gex))

pain_scores_low <- pain_scores[rownames(gex),]

df <- data.frame(x = gex$TF,
                 y = pain_scores_low$pain_scores,
                 ps = pain_scores_low$pain_scores)

low.kendall.result <- cor.test(df$x, df$y, method="kendall")
low.pearson.result <- cor.test(df$x, df$y, method="pearson")
low.spearman.result <- cor.test(df$x, df$y, method="spearman")

tf_vs_pain<-ggplot(df, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="orange", size=7) +
  #geom_point(aes(size=10), show.legend=FALSE,colour="black") +
  geom_smooth(method = "lm", se = FALSE, color="black") +
  theme(#text = element_text(size = 30),
    #axis.text.x = element_text(size = 30),
    axis.title.x = element_blank(),
    #axis.text.y = element_text(size = 30),
    axis.title.y = element_blank())+
  labs(x=NULL, y=NULL)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))



print("Low inflammation patient HOOS/KOOS result")

print(paste('Kendall Correlation test P-value:',low.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',low.kendall.result$estimate))

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_tf_expression_vs_pain_scores_lowinf",".png",sep=""),
       width = 7, height = 10, dpi = 300, units = "in", device='png')

#==========================================================================================#
#==============================================================================================================================#

# Low inflammatory patients

gex <- read.csv(paste(input_dir, 'GEX_TF_bulk_highinf.csv',sep=""), row.names = 1)

gex <- data.frame(t(gex))

pain_scores_high <- pain_scores[rownames(gex),]

df <- data.frame(x = gex$TF,
                 y = pain_scores_high$pain_scores,
                 ps = pain_scores_high$pain_scores)

high.kendall.result <- cor.test(df$x, df$y, method="kendall")
high.pearson.result <- cor.test(df$x, df$y, method="pearson")
high.spearman.result <- cor.test(df$x, df$y, method="spearman")

tf_vs_pain<-ggplot(df, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="gray", size=7) +
  #geom_point(aes(size=10), show.legend=FALSE,colour="black") +
  geom_smooth(method = "lm", se = FALSE, color="black") +
  theme(#text = element_text(size = 30),
    #axis.text.x = element_text(size = 30),
    axis.title.x = element_blank(),
    #axis.text.y = element_text(size = 30),
    axis.title.y = element_blank())+
  labs(x=NULL, y=NULL)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_text(size=25))



print("High inflammation patient HOOS/KOOS result")

print(paste('Kendall Correlation test P-value:',high.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',high.kendall.result$estimate))

ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_tf_expression_vs_pain_scores_highinf",".png",sep=""),
       width = 7, height = 10, dpi = 300, units = "in", device='png')