# Author: Komal S. Rathi
# Date: 04/28/2020
# Function: Calculate distances between strains

library(usedist)
library(edgeR)
library(tidyverse)
library(xlsx)
library(reshape2)

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
  
  # if group is NULL, do not filter
  if(!is.null(group)){
    metaDat <- metaDat %>%
      filter(label2 == group) %>%
      column_to_rownames("sample")
  } else {
    metaDat <- metaDat %>%
      mutate(strain = paste0(strain, '-', label2)) %>%
      column_to_rownames("sample")
  }
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

# all mice together (4 strains x 2 groups)
res <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = NULL)
res <- melt(res)
res <- res %>%
  filter(value != 0)

# calculate distances between strains - non-exercised mice
# B is close to E < I < A
# E is close to I < B < A
# I is close to E < B < A
nonex <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = "non_exercised")
write.xlsx(nonex,  file = 'results/qc/strain-distances.xlsx', sheetName = 'nonex', append = TRUE)
# calculate distances between strains - exercised mice
# B is close to I < E < A
# E is close to I < B < A
# I is close to E < B < A
ex <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = "exercised")
write.xlsx(ex,  file = 'results/qc/strain-distances.xlsx', sheetName = 'ex', append = TRUE)

# so the only change we see is the B strains were close to E strains in nonexercised
# but become closer to I strain in exercised mice

# now we will do pairwise comparison between nonex and ex
nonex <- melt(nonex)
nonex <- nonex %>%
  filter(value != 0) %>%
  mutate(label = "nonex")
ex <- melt(ex)
ex <- ex %>%
  filter(value != 0) %>%
  mutate(label = "ex")
total <- rbind(ex, nonex)

# 1. check if data is paired
# 2. sample size < 30, so let's do Shapiro Test
d <- with(total, 
          value[label == "ex"] - value[label == "nonex"])
shapiro.test(d) # Shapiro-Wilk normality test for the differences
# p-value > 0.05 
# implying that the distribution of the differences (d) are not significantly different from normal distribution. 
# In other words, we can assume the normality.
# 3. So, now we can perform paired t-test
# if the p-value was significant, we would have to do KS test
total2 <- dcast(total, Var1+Var2~label, value.var = 'value')
res <- t.test(value~label, data = total, paired = TRUE) # same
res2 <- t.test(total2$nonex, total2$ex, paired = TRUE) # same

# Conclusion
# Pvalue > 0.05. 
# We can then reject null hypothesis and conclude that the distances between strains for nonexercised is not significantly different from exercised mice.


