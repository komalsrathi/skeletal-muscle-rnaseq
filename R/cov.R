library(usedist)
library(edgeR)
library(tidyverse)
library(xlsx)
library(reshape2)
library(bioplotr)

source('R/filterExpr.R')

load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$litter[meta$litter == ""] <- "U"
rownames(meta) <- meta$sample
meta$label2 <- ifelse(meta$label == "non_exercised", meta$label, "exercised")
expr.counts.mat <- expr.counts.mat[,rownames(meta)]

# filter expression
expr.counts.mat <- filterExpr(expr.counts.mat)

# create design
var <- factor(meta[,'strain'])
design <- model.matrix(~0+var)
colnames(design) <- levels(var)
rownames(design) <- meta$sample

# voom normalize data
y <- DGEList(counts = as.matrix(expr.counts.mat), genes = rownames(expr.counts.mat))
y <- calcNormFactors(y)
v <- voom(counts = y, design = design, plot = FALSE)
voomData <- v$E

# covariates
covariates <- read.xlsx('data/Covariable list Patrick Schaefer - RNASeq soleus.xls', sheetIndex = 3)
covariates <-  covariates %>%
  column_to_rownames("mice")
covariates <- covariates[rownames(meta),]
covariates[is.na(covariates)] <- 0
covariates <- covariates %>%
  as_tibble()

# plot 
plot_drivers(dat = voomData, clin = covariates, block = "VO2max", unblock = "running.time", p_adj = "BH")

