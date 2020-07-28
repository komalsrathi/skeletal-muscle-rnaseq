# Author: Komal S. Rathi
# Date: 03/05/2020
# Function: Clustering t-SNE
# combat adjust - not needed because library prep was done together

library(optparse)
library(tidyverse)
library(reshape2)
library(Rtsne)
library(sva)
library(ggpubr)
source('~/Projects/Utils/pubTheme.R')

option_list <- list(
  make_option(c("--fpkm_matrix"), type = "character",
              help = "RData object of FPKM"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file for samples (.tsv)"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
fpkm_matrix <- opt$fpkm_matrix
annot <- opt$annot
outdir <- opt$outdir
meta_file <- opt$meta_file

# load meta
df <- read.delim(meta_file, stringsAsFactors = F)

# load expression and remove all zeros
load(fpkm_matrix)
expr.fpkm <- expr.fpkm[[1]]
expr.fpkm <- expr.fpkm[apply(expr.fpkm[,-1], 1, function(x) !all(x==0)),]

# meta
df$type <- gsub("[0-9].*", "", df$sample)
rownames(df) <- df$sample
df <- df[colnames(expr.fpkm),]

# house keeping genes
hkgenes <- c("Actb", "Tuba1a", "Gapdh", "Ldha", "Rpl19")

# correct for strain
expr.fpkm.corrected <- ComBat(dat = log2(expr.fpkm + 1), batch = df$type) 

# t-SNE before combat adjustment
set.seed(42)
tsneOut <- Rtsne(t(log2(expr.fpkm + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
p <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(color = label)) +
  geom_text(aes(label = type)) + 
  theme_bw() +
  ggtitle("T-SNE Clustering (Before Batch Correction)") +
  theme_Publication2() + xlab("PC1") + ylab("PC2") 

expr.fpkm.hk <- expr.fpkm[rownames(expr.fpkm) %in% hkgenes,]
tsneOut <- Rtsne(t(log2(expr.fpkm.hk + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
q <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(color = label)) +
  geom_text(aes(label = type)) + 
  theme_bw() +
  ggtitle("T-SNE Clustering HK genes (Before Batch Correction)") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

# t-SNE after combat adjustment
set.seed(42)
tsneOut <- Rtsne(t(log2(expr.fpkm.corrected + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
r <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(color = label)) +
  geom_text(aes(label = type)) + 
  theme_bw() +
  ggtitle("T-SNE Clustering (After Batch Correction)") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

expr.fpkm.corrected.hk <- expr.fpkm.corrected[rownames(expr.fpkm.corrected) %in% hkgenes,]
tsneOut <- Rtsne(t(log2(expr.fpkm.corrected.hk + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
s <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(color = label)) +
  geom_text(aes(label = type)) + 
  theme_bw() +
  ggtitle("T-SNE Clustering HK genes (After Batch Correction)") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

fname <- file.path(outdir, 'tSNE_plot.pdf')
ggarrange(p, q, r, s, common.legend = T) %>%
  ggexport(filename = fname, width = 12, height = 10)

