# Author: Komal S. Rathi
# Function: Batch effect correction using atrial cardiomyocytes expression
# Date: 08/19/2020

library(tidyverse)
library(sva)

# Expression plot of heart-specific genes
genes <- read.delim('data/heart/heart-specific-genes.txt', stringsAsFactors = F)
genes <- unique(genes)
gene.sym <- unique(genes$gene_symbol)

# FPKM expression
load('data/heart/heart_collapsed_fpkm_matrix.RData')
expr.fpkm <- expr.fpkm[[1]]
expr.fpkm <- expr.fpkm[rownames(expr.fpkm) %in% gene.sym,]
expr.fpkm <- melt(as.matrix(expr.fpkm), varnames = c("gene_symbol", "sample"), value.name = 'fpkm')

# meta data
meta <- read.delim('data/heart/heart-meta-data.txt', stringsAsFactors = F)
expr.fpkm <- expr.fpkm %>%
  inner_join(meta, by = c('sample'))

# merge grouping of genes
expr.fpkm <- expr.fpkm %>%
  inner_join(genes, by = c('gene_symbol')) %>%
  filter(cell_type != "")  

# only atrial cardiomyocytes
celltypes <- unique(expr.fpkm$cell_type)
dat <- expr.fpkm %>%
  filter(cell_type %in% celltypes[1]) %>%
  mutate(fpkm = log2(fpkm + 1))

# define high and low batches
batch  <- dat %>%
  group_by(sample) %>%
  summarise(median = median(fpkm)) %>%
  arrange(median)
batch <- batch %>%
  mutate(batch = ifelse(median >= mean(median), 'High', 'Low'))

# add batch variable to metadata
meta <- meta %>%
  inner_join(batch, by = 'sample') %>%
  dplyr::select(-c(median)) %>%
  mutate(tmp = sample) %>%
  column_to_rownames('tmp')

# load counts
load('data/heart/heart_collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]

# now batch correct and back transform
expr.counts.mat <- expr.counts.mat[,rownames(meta)]
identical(colnames(expr.counts.mat), rownames(meta))
var <- factor(meta[,'strain'])
design <- model.matrix(~0+var)
colnames(design) <- levels(var)
rownames(design) <- meta$sample
batch <- factor(meta[,'batch'])
expr.counts.mat.bc <- ComBat(log2(expr.counts.mat + 1), batch = batch, par.prior = T)

# back transform
expr.counts.mat.bc <- 2^expr.counts.mat.bc
expr.counts <- list(expr.counts.mat.bc, expr.counts.annot)
save(expr.counts, file = 'data/heart/heart_collapsed_counts_matrix_batchcorrected.RData')

