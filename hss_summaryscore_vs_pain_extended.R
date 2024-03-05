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

pain_scores <- read.csv(paste(input_dir, 'pain_scores_hooskoos.csv', sep=""), row.names=1)

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_pain815.csv',sep=""), row.names = 1)

patients_hooskoos <- intersect(row.names(pain_scores), colnames(gex_zscores))

pain_scores_hooskoos <- pain_scores[patients_hooskoos,]
gex_zscores <- gex_zscores[, patients_hooskoos]

set.seed(42)
gex_zscores.umap = umap(t(gex_zscores))

df <- data.frame(x = gex_zscores.umap$layout[,1],
                 y = gex_zscores.umap$layout[,2],
                 ps = pain_scores_hooskoos$pain_scores)

sp2<-ggplot(df, aes(x, y, color=ps)) +
  geom_point(aes(size=ps))
sp2+scale_color_gradient(low="red", high="green")

set.seed(42)
perp=11 # Largest possible rounded degree according to pain-score graph. Graph degrees: 11.762489490419542
gex_zscores.tsne <- Rtsne(t(gex_zscores), dims = 1, perplexity=perp, verbose=F) # NOTE: perplexity needs to be consistent!

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                              y = pain_scores_hooskoos$pain_scores,
                              ps = pain_scores_hooskoos$pain_scores)

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

all.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
all.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
all.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

df_sampled.tsne <- data.frame(x = -sign(as.numeric(kendall.result$estimate[1]))*gex_zscores.tsne$Y,
                              y = pain_scores$pain_scores,
                              ps = pain_scores$pain_scores)

all.kendall.result <- cor.test(df_sampled.tsne$x, df.tsne$y, method="kendall")
all.pearson.result <- cor.test(df_sampled.tsne$x, df.tsne$y, method="pearson")
all.spearman.result <- cor.test(df_sampled.tsne$x, df.tsne$y, method="spearman")

sp_ls.tsne<-ggplot(df_sampled.tsne, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="blue", size=7) +
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



print("All-patient HOOS/KOOS result")
print(paste("Perplexity:", perp))
print(paste('Kendall Correlation test P-value:',all.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',all.kendall.result$estimate))

figure_title <- paste('Kendall Correlation test P-value:',format(kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne + 
  scale_y_continuous(limits=c(0, 60)) +  #name="HOOS/KOOS Pain Score", 
  labs(x=NULL, y=NULL)
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_summary_scores_vs_pain_scores_highlowinfla_pain_genes_paper_Fig2_p11",".png",sep=""),
       width = 7, height = 10, dpi = 150, units = "in", device='png')
#dev.off()

#######################################################################################

set.seed(42)
perp = 7 # Following low inflammation
gex_zscores.tsne <- Rtsne(t(gex_zscores), dims = 1, perplexity=perp, verbose=F) # NOTE: perplexity needs to be consistent!

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                      y = pain_scores_hooskoos$pain_scores,
                      ps = pain_scores_hooskoos$pain_scores)

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

all.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
all.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
all.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

df_sampled.tsne <- data.frame(x = -sign(as.numeric(kendall.result$estimate[1]))*gex_zscores.tsne$Y,
                              y = pain_scores$pain_scores,
                              ps = pain_scores$pain_scores)

all.kendall.result <- cor.test(df_sampled.tsne$x, df.tsne$y, method="kendall")
all.pearson.result <- cor.test(df_sampled.tsne$x, df.tsne$y, method="pearson")
all.spearman.result <- cor.test(df_sampled.tsne$x, df.tsne$y, method="spearman")

sp_ls.tsne<-ggplot(df_sampled.tsne, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="blue", size=7) +
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



print("All-patient HOOS/KOOS result")
print(paste("Perplexity:", perp))
print(paste('Kendall Correlation test P-value:',all.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',all.kendall.result$estimate))

figure_title <- paste('Kendall Correlation test P-value:',format(kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne + 
  scale_y_continuous(limits=c(0, 60)) +  #name="HOOS/KOOS Pain Score", 
  labs(x=NULL, y=NULL)
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_summary_scores_vs_pain_scores_highlowinfla_pain_genes_paper_Fig2_p7",".png",sep=""),
       width = 7, height = 10, dpi = 150, units = "in", device='png')
#dev.off()


#######################################################################################
# HOOS/KOOS High Inflammation Patients
#######################################################################################

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_pain815_highinf.csv',sep=""), row.names = 1)

patients_high <- colnames(gex_zscores)

pain_scores_highinf <- pain_scores[patients_high,]
gex_zscores <- gex_zscores[, patients_high]

set.seed(42)
perp=3 # Largest possible rounded degree according to pain-score graph. Graph degrees: 4.840544407564654
gex_zscores.tsne <- Rtsne(t(gex_zscores), dims = 1, perplexity=perp, verbose=F) # NOTE: perplexity needs to be consistent!

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                      y = pain_scores_highinf$pain_scores,
                      ps = pain_scores_highinf$pain_scores)

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

highinf.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
highinf.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
highinf.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

# print("High inflammation-patient HOOS/KOOS result")
# print(paste("Perplexity:", perp))
# print(paste('Kendall Correlation test P-value:',highinf.kendall.result$p.value))
# print(paste('Kendall Correlation coefficient:',highinf.kendall.result$estimate))

gex_zscores.tsne$Y <- -sign(highinf.kendall.result$estimate)*gex_zscores.tsne$Y

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                      y = pain_scores_highinf$pain_scores,
                      ps = pain_scores_highinf$pain_scores)

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

highinf.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
highinf.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
highinf.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

sp_ls.tsne<-ggplot(df.tsne, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="green", size=7) +
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


print("High inflammation-patient HOOS/KOOS result")
print(paste("Perplexity:", perp))
print(paste('Kendall Correlation test P-value:',highinf.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',highinf.kendall.result$estimate))

figure_title <- paste('Kendall Correlation test P-value:',format(highinf.kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne + 
  scale_y_continuous(limits=c(0, 60)) +  #name="HOOS/KOOS Pain Score", 
  labs(x=NULL, y=NULL)
ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_summary_scores_vs_pain_scores_highinfla_pain_genes_paper_Fig2_p3",".png",sep=""),
       width = 7, height = 10, dpi = 150, units = "in", device='png')
#dev.off()

########################

gex_zscores <- read.csv(paste(input_dir, 'GEX_zscores_bulk_pain815_lowinf.csv',sep=""), row.names = 1)

patients_low <- colnames(gex_zscores)

pain_scores_lowinf <- pain_scores[patients_low,]
gex_zscores <- gex_zscores[, patients_low]

set.seed(42)
perp=7 # 
gex_zscores.tsne <- Rtsne(t(gex_zscores), dims = 1, perplexity=perp, verbose=F) # NOTE: perplexity needs to be consistent!

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                      y = pain_scores_lowinf$pain_scores,
                      ps = pain_scores_lowinf$pain_scores)

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

lowinf.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
lowinf.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
lowinf.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

# print("High inflammation-patient HOOS/KOOS result")
# print(paste("Perplexity:", perp))
# print(paste('Kendall Correlation test P-value:',highinf.kendall.result$p.value))
# print(paste('Kendall Correlation coefficient:',highinf.kendall.result$estimate))

gex_zscores.tsne$Y <- -sign(lowinf.kendall.result$estimate)*gex_zscores.tsne$Y

df.tsne <- data.frame(x = gex_zscores.tsne$Y,
                      y = pain_scores_lowinf$pain_scores,
                      ps = pain_scores_lowinf$pain_scores)

kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")

lowinf.kendall.result <- cor.test(df.tsne$x, df.tsne$y, method="kendall")
lowinf.pearson.result <- cor.test(df.tsne$x, df.tsne$y, method="pearson")
lowinf.spearman.result <- cor.test(df.tsne$x, df.tsne$y, method="spearman")

print("Low inflammation-patient HOOS/KOOS result")
print(paste("Perplexity:", perp))
print(paste('Kendall Correlation test P-value:',lowinf.kendall.result$p.value))
print(paste('Kendall Correlation coefficient:',lowinf.kendall.result$estimate))

sp_ls.tsne<-ggplot(df.tsne, aes(x, y, label=ps)) +
  geom_point(aes(size=15), show.legend=FALSE,colour="blue", size=7) +
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



figure_title <- paste('Kendall Correlation test P-value:',format(highinf.kendall.result$p.value,scientific=TRUE))
#+ggtitle(figure_title)
sp_ls.tsne + 
  scale_y_continuous(limits=c(0, 60)) +  #name="HOOS/KOOS Pain Score", 
  labs(x=NULL, y=NULL)
#ggsave(paste("/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/hss_summary_scores_vs_pain_scores_low_infla_genes_paper_Fig2",".png",sep=""))
#dev.off()



