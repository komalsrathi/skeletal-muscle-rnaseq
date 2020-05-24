# Author: Komal S. Rathi
# Date: 04/28/2020
# Function: Calculate distances between strains

library(usedist)
library(edgeR)
library(tidyverse)
library(xlsx)
library(reshape2)

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

calc.dist <- function(exprDat, metaDat, group, groupSamples = TRUE) {
  
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
  
  if(groupSamples == TRUE){
    mat.dist <- as.matrix(dist_multi_centroids(dm4, g4))
  } else {
    mat.dist <- dist_groups(dm4, g4)
  }
  return(mat.dist)
}

# 1. Samples Group by Strain 
test.fun <- function(exprDat, metaDat, fname){
  # non-exercised/exercised
  nonex <- calc.dist(exprDat, metaDat, group = "non_exercised")
  ex <- calc.dist(exprDat, metaDat, group = "exercised")
  if(!file.exists(fname)){
    write.xlsx(nonex, file = fname, sheetName = "nonex", append = T)
    write.xlsx(ex, file = fname, sheetName = "ex", append = T)
  }
  
  # melt matrix
  # non-exercised
  nonex <- melt(nonex)
  nonex <- nonex %>%
    filter(value != 0) %>%
    mutate(label = "nonex")
  
  # exercised
  ex <- melt(ex)
  ex <- ex %>%
    filter(value != 0) %>%
    mutate(label = "ex")
  total <- rbind(ex, nonex)
  
  # shapiro test  
  d <- with(total, value[label == "ex"] - value[label == "nonex"])
  st <- shapiro.test(d)$p.value # Shapiro-Wilk normality test for the differences
  print(paste0("Test of Normality Pval:", round(st, 2)))
  if(st < 0.05){
    print("Wilcoxon Test")
    res <-  wilcox.test(value~label, data = total, paired = TRUE, correct=FALSE)$p.value
  } else {
    print("Paired T-test")
    res <- t.test(value~label, data = total, paired = TRUE)$p.value # same
  }
  print(res)
}
# test.fun(exprDat = expr.counts.mat, metaDat = meta, fname = "results/distances/strain-distances-Counts.xlsx") # Counts 
# Pvalue > 0.05. Distances between strains for non-exercised are not significantly different from exercised mice.

test.fun(exprDat = voomData, metaDat = meta, fname  = "results/distances/distances-across-strains.xlsx")  # Voom normalized
# Pvalue < 0.05. Distances between strains for non-exercised are significantly different from exercised mice.

# 2. Between Sample Comparison
nonex <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = "non_exercised", groupSamples = FALSE)
ex <- calc.dist(exprDat = expr.counts.mat, metaDat = meta, group = "exercised", groupSamples = FALSE)
nonex <- nonex %>%
  as.data.frame() %>%
  mutate(Label2 = "nonex")
ex <- ex %>%
  as.data.frame() %>%
  mutate(Label2 = "ex")
total <- rbind(ex, nonex)
total <- total[grep('Within', total$Label),]
total <- total[,c('Ex-Label','Ex-Distance','Nonex-Distance')]
return.p <- function(x) {
  print(x)
  tt <- t.test(Distance~Label2, data = x, paired = FALSE)$p.value
  wt <- wilcox.test(Distance~Label2, data = x, paired = FALSE)$p.value
  res <- data.frame(t.test = tt, wilcox.test = wt)
  return(res)
}
res <- plyr::ddply(total, 
            .variables = 'Label', 
            .fun = function(x) return.p(x))
write.xlsx(x = res, file = 'results/distances/distances-between-strains.xlsx', row.names = FALSE)
            