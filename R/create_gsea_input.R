# Author: Komal S. Rathi
# Function: Create GSEA input files
# Date: 05/26/2020

# libraries
library(tidyverse)
library(reshape2)

# load count matrix
load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.mat <- expr.counts.mat[apply(expr.counts.mat!=0, 1, all),]

#expr.counts.annot <-  expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$litter[meta$litter == ""] <- "U"
rownames(meta) <- meta$sample
meta$label2 <- ifelse(meta$label == "non_exercised", meta$label, "exercised")

# function
# 1. Overall Ex vs Nonex: So only 2 groups
# 2. By strain Ex vs Nonex: So 8 groups
system('mkdir -p results/gsea')

gsea.input <- function(counts.collapsed, meta, groupby, gct.file, cls.file) {
  # make file GSEA web version
  # expression file  -  counts file 
  # df <- data.frame(snames = colnames(counts.collapsed), colsplit(colnames(counts.collapsed), pattern = "_", names = c("sample","type")))
  if(groupby == "strain"){
    meta$strain <- gsub(' ','_',meta$strain)
    meta$label2 <- paste0(meta$strain, '_',meta$label2)
  }
  
  # order
  meta <- meta[order(meta$label2),]
  counts.collapsed <- counts.collapsed[,rownames(meta)]
  
  gct <- counts.collapsed[,rownames(meta)]
  add <- data.frame(NAME = c("#1.2",nrow(gct),"NAME"), Description = c('', ncol(gct), "Description"))
  total.cols <- ncol(gct) + 2
  add[,3:total.cols] <- ''
  colnames(add)[3:total.cols] <- colnames(gct)
  add[3,3:total.cols] <- colnames(gct)
  annot <- data.frame(NAME = rownames(gct), Description = 'na')
  annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
  add <- rbind(add, annot)
  write.table(add, file = gct.file, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # phenotype file
  groups <- levels(factor(meta$label2))
  ngroups <- length(groups)
  ph <- matrix(nrow = 3, ncol = ncol(gct))
  # first row
  ph[1,1] <- ncol(gct)
  ph[1,2] <- ngroups
  ph[1,3] <- 1
  # second row
  ph[2,1:ngroups] <- groups
  ph[2,1] <- paste0('# ', ph[2,1])
  # third row
  ph[3,] <- meta$label2
  ph <- as.data.frame(ph)
  write.table(ph, file = cls.file, quote = F, sep = " ", na = "", col.names = F, row.names = F)
  
  # chp file
  # chp <- data.frame('Probe Set ID' = rownames(gct), 'Gene Symbol' = rownames(gct), 'Gene Title' = 'na', check.names = F)
  # write.table(chp, file = 'data/gsea_input/mapping.chip', quote = F, sep = "\t", row.names = F)
}

# by group
gsea.input(counts.collapsed = expr.counts.mat,
           meta = meta,
           groupby = 'label2',
           gct.file = 'results/gsea/matrix_byGroup.gct',
           cls.file = 'results/gsea/phenotype_byGroup.cls')

# by strain
gsea.input(counts.collapsed = expr.counts.mat,
           meta = meta,
           groupby = 'strain',
           gct.file = 'results/gsea/matrix_byStrain.gct',
           cls.file = 'results/gsea/phenotype_byStrain.cls')
