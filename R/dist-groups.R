# Author: Komal S. Rathi
# Date: 04/28/2020
# Function: Calculate distances between strains

library(usedist)
library(edgeR)
library(tidyverse)
library(xlsx)

m4 <- matrix(1:16, nrow=4, dimnames=list(LETTERS[1:4]))
dm4 <- dist(m4)
g4 <- rep(c("Control", "Treatment"), each=2)
dist_groups(dm4, g4)
dist_multi_centroids(dm4, g4)
dist_to_centroids(dm4, g4)


load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$litter[meta$litter == ""] <- "U"
rownames(meta) <- meta$sample
meta$label2 <- ifelse(meta$label == "non_exercised", meta$label, "exercised")
expr.counts.mat <- expr.counts.mat[,rownames(meta)]


calc.dist <- function(exprDat, metaDat, group) {
  metaDat <- metaDat %>%
    filter(label2 == group) %>%
    column_to_rownames("sample")
  exprDat <- exprDat[,rownames(metaDat)]

  # filter low expression
  keep.exprs <- filterByExpr(exprDat)
  exprDat <- exprDat[keep.exprs,]
  print(dim(exprDat))

  # calculate distance
  dm4 <- dist(t(exprDat))
  g4 <- metaDat$strain

  mat.dist <- as.matrix(dist_multi_centroids(dm4, g4))
  return(mat.dist)
}

# calculate distances between strains - non-exercised mice
# B is close to E < I < A
# E is close to I < B < A
# I is close to E < B < A
res <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = "non_exercised")
write.xlsx(res,  file = 'results/qc/strain-distances.xlsx', sheetName = 'nonex', append = TRUE)
# calculate distances between strains - exercised mice
# B is close to I < E < A
# E is close to I < B < A
# I is close to E < B < A
res <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = "exercised")
write.xlsx(res,  file = 'results/qc/strain-distances.xlsx', sheetName = 'ex', append = TRUE)

# so the only change we see is the B strains were close to E strains in nonexercised
# but become closer to I strain in exercised mice
