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

library(gprofiler2)

library(openxlsx)

data_dir <- '/Users/baizilong/Documents/Dana-dataset/RA_paper_updated/'

df_nonpain_genes <- read.csv(paste(data_dir, 'nonpain_low_inflammation_gene_symbols.csv',sep=""), header=FALSE)
nonpain_genes <- df_nonpain_genes$V1

df_background_genes <- read.csv(paste(data_dir, 'bulk_background_gene_symbols.csv',sep=""), header=FALSE)
background_genes <- df_background_genes$V1

gostres <- gost(nonpain_genes, significant = FALSE, correction_method = "gSCS", domain_scope = "custom", custom_bg = background_genes)

# substring = "neuro"
# gostres$result$term_name
# 
# matching_elements <- grep(substring, gostres$result$term_name, value = TRUE, ignore.case = TRUE)
# matching_elements

PEA_full_results <- data.frame(gostres$result)
write.xlsx(PEA_full_results, file = paste(data_dir,"nonpain_PEA_full.xlsx",sep=""), row.names = FALSE)
# interest_pathways <- c("Neurogenesis", "Nervous system process", "Nervous system development", "Multicellular organism development", "Generation of neurons", 
#                        "Developmental process", "Cell morphogenesis involved in neuron differentiation", "Anatomical structure morphogenesis",
#                        "Anatomical structure development")
# nonpain_PEA <- PEA_results[tolower(PEA_results$term_name) %in% tolower(interest_pathways), ]
# dim(nonpain_PEA)
# 
# nonpain_PEA[c('p_value', 'term_name', 'term_id')]


#====
  
interest_pathways <- c("neurogenesis", "nervous system process", "nervous system development", "multicellular organism development", "generation of neurons", 
                         "developmental process", "cell morphogenesis involved in neuron differentiation", "anatomical structure morphogenesis",
                         "anatomical structure development")
nonpain_PEA_selected <- as.data.frame(PEA_full_results[PEA_full_results$term_name %in% interest_pathways, ])
nonpain_PEA_selected <- data.frame(nonpain_PEA_selected)
dim(nonpain_PEA_selected)
rownames(nonpain_PEA_selected) <- nonpain_PEA_selected$term_name

write.xlsx(nonpain_PEA_selected, file = paste(data_dir,"nonpain_PEA_selected.xlsx",sep=""), row.names = FALSE)

# nonpain_PEA[c('p_value', 'term_name', 'term_id')]
#nonpain_PEA_save <- apply(nonpain_PEA,2,as.character) #https://stackoverflow.com/questions/24829027/unimplemented-type-list-when-trying-to-write-table
#write.csv(nonpain_PEA_save, file=paste(data_dir,"nonpain_PEA_selected.csv",sep=""), row.names = FALSE)

#========= READ pain PEA
df_pain_PEA <- read.csv(paste(data_dir, 'Supplementary Table 1 Pathway Analysis.csv',sep="")) 
pain_PEA <- df_pain_PEA[df_pain_PEA$term_name %in% interest_pathways,]

rownames(pain_PEA) <- pain_PEA$term_name

pain_PEA <- pain_PEA[interest_pathways,]
nonpain_PEA_selected <- nonpain_PEA_selected[interest_pathways,]

nonpain_PEA_vals <- -log10(nonpain_PEA_selected$p_value)
pain_PEA_vals <- pain_PEA$negative_log10_of_adjusted_p_value
pathways <- c("Neurogenesis", "Nervous system process", "Nervous system development", "Multicellular organism development", "Generation of neurons", 
              "Developmental process", "Cell morphogenesis involved in neuron differentiation", "Anatomical structure morphogenesis",
              "Anatomical structure development")

painVSnonpain <- data.frame(pathways = pathways, Nonpain = nonpain_PEA_vals, Pain = pain_PEA_vals)

write.csv(painVSnonpain, file=paste(data_dir,"painVSnonpain_PEA.csv",sep=""), row.names = FALSE)
