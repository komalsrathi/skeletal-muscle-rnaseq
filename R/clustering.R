# Author: Komal S. Rathi
# Date: 03/05/2020
# Function: Clustering using UMAP, PCA and t-SNE

library(optparse)
library(tidyverse)
library(reshape2)
library(Rtsne)
library(sva)
library(ggpubr)
library(ggplot2)
library(matrixStats)
library(umap)
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
outdir <- opt$outdir
meta_file <- opt$meta_file

# make output directories
outdir <- file.path(outdir, 'clustering-qc')
dir.create(outdir, showWarnings = F, recursive = TRUE)

# load meta
df <- read.delim(meta_file, stringsAsFactors = F)

# load expression and remove all zeros
load(fpkm_matrix)
expr.fpkm <- expr.fpkm[[1]]
expr.fpkm <- expr.fpkm[apply(expr.fpkm[,-1], 1, function(x) !all(x==0)),]

# meta
rownames(df) <- df$sample
df <- df[colnames(expr.fpkm),]

# house keeping genes
hkgenes <- read.delim('data/gene-lists/housekeeping-genes-mouse.txt', stringsAsFactors = F)
hkgenes <- hkgenes$genes

# define input files
# house keeping genes
expr.fpkm.hk <- expr.fpkm[rownames(expr.fpkm) %in% hkgenes,]
# 500 most variable genes
rv <- rowVars(as.matrix(expr.fpkm))
select <- order(rv, decreasing=TRUE)[seq_len(500)]
expr.fpkm.most.var <- expr.fpkm[select,]

# t-SNE 
set.seed(42)
tsneOut <- Rtsne(t(log2(expr.fpkm + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
p <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(shape = label, color = strain)) +
  theme_bw() +
  ggtitle("T-SNE Clustering") +
  theme_Publication2() + xlab("PC1") + ylab("PC2") 

set.seed(42)
tsneOut <- Rtsne(t(log2(expr.fpkm.hk + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
q <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(shape = label, color = strain)) +
  theme_bw() +
  ggtitle("T-SNE Clustering HK genes") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

set.seed(42)
tsneOut <- Rtsne(t(log2(expr.fpkm.most.var + 1)), initial_dims = 50, perplexity = 8, max_iter = 1000)
tsneOut <- data.frame(tsneOut$Y, df)
r <- ggplot(tsneOut, aes(X1, X2)) +
  geom_point(size = 5, alpha = 0.5, aes(shape = label, color = strain)) +
  theme_bw() +
  ggtitle("T-SNE Clustering Top 500 Most Variable Genes") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

fname <- file.path(outdir, 'tSNE_plot.pdf')
ggarrange(p, q, r, common.legend = T) %>%
  ggexport(filename = fname, width = 20, height = 10)

# PCA
prData <- prcomp(log2(expr.fpkm + 1))
pca.data <- prData$rotation
pca.data <- data.frame(pca.data)[1:4]
pca.data <- data.frame(pca.data, df)
p <- ggplot(pca.data, aes(PC1, PC2)) +
  geom_point(size = 5, alpha = 0.5, aes(shape = label, color = strain)) +
  theme_bw() +
  ggtitle("PCA Clustering") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

prData <- prcomp(log2(expr.fpkm.hk + 1))
pca.data <- prData$rotation
pca.data <- data.frame(pca.data)[1:4]
pca.data <- data.frame(pca.data, df)
q <- ggplot(pca.data, aes(PC1, PC2)) +
  geom_point(size = 5, alpha = 0.5, aes(shape = label, color = strain)) +
  theme_bw() +
  ggtitle("PCA Clustering HK genes") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

prData <- prcomp(log2(expr.fpkm.most.var + 1))
pca.data <- prData$rotation
pca.data <- data.frame(pca.data)[1:4]
pca.data <- data.frame(pca.data, df)
r <- ggplot(pca.data, aes(PC1, PC2)) +
  geom_point(size = 5, alpha = 0.5, aes(shape = label, color = strain)) +
  theme_bw() +
  ggtitle("PCA Clustering Top 500 Most Variable Genes") +
  theme_Publication2() + xlab("PC1") + ylab("PC2")

fname <- file.path(outdir, 'pca_plot.pdf')
ggarrange(p, q, r, common.legend = T) %>%
  ggexport(filename = fname, width = 20, height = 10)

# pheatmap
library(pheatmap)
#pdf(file = file.path(outdir, "heatmap.pdf"), width = 10, height = 15)
pheatmap(log2(expr.fpkm.most.var + 1), 
         scale = "row", show_rownames = F, 
         main = "Top 500 most variable genes",
         annotation_col = df[,c("strain", "label")], 
         filename = file.path(outdir, "heatmap.pdf"), width = 15, height = 10)
#dev.off()
