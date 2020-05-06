# Author: Komal S. Rathi
# Date: 03/05/2020
# Function: Clustering t-SNE
# combat adjust - not needed because library prep was done together

setwd('~/Projects/Skeletal_Muscle_Patrick/')

library(tidyverse)
library(reshape2)
library(Rtsne)
source('~/Projects/Utils/pubTheme.R')

load('data/skeletal_muscle_mm10.RData')
annot <- data.table::fread('data/gencode.vM17.annotation.txt')
annot <- annot %>%
  select(c(gene_id, gene_symbol, biotype)) %>%
  unique()

# remove all zeros
expr.fpkm <- expr.fpkm[apply(expr.fpkm[,-1], 1, function(x) !all(x==0)),]
rownames(expr.fpkm) <- expr.fpkm$gene_id
expr.fpkm$gene_id <- NULL

# meta
df <- read.delim('data/meta-data.txt', stringsAsFactors = F)
df$type <- gsub("[0-9].*", "", df$sample)
rownames(df) <- df$sample
df <- df[colnames(expr.fpkm),]

# t-SNE before combat adjustment
set.seed(42)
tsneOut <- Rtsne(t(expr.fpkm), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
p <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(color = label)) +
  geom_text(aes(label = type)) + 
  theme_bw() +
  ggtitle("T-SNE Clustering (Before Batch Correction)") +
  theme_Publication2() + xlab("PC1") + ylab("PC2") 
p
ggsave(p, filename = 'results/tSNE_plot.pdf', width = 8, height = 6)

