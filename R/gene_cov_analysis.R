# Author: Komal S. Rathi
# Function: Gene expression correlation with covariates analysis
# 1. Overall nonex vs ex
# 2. Per Strain; nonex vs ex
library(usedist)
library(edgeR)
library(tidyverse)
library(xlsx)
library(reshape2)
library(bioplotr)
library(broom)
library(ggpubr)

source('R/filterExpr.R')
source('../Utils/pubTheme.R')

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

# merge covariates to metadata
if(identical(rownames(covariates), rownames(meta))){
  meta <- cbind(meta, covariates)
}

# transpose matrix
voomData.trn <- t(voomData)
if(identical(rownames(voomData.trn), rownames(meta))){
  tmp <- cbind(meta, voomData.trn[,1:2])
  tmp <- tmp %>%
    dplyr::select(c(label2, strain, VO2max, `mt-Co1`)) %>%
    mutate(gene = `mt-Co1`)
}

# boxplots of covariates
for(i in 6:ncol(meta)){
  forbox <- meta
  ylab <- colnames(forbox)[i]
  colnames(forbox)[i]  <- 'covariate'
  fname <- paste0('results/covars/',ylab,'.png')
  my_comparisons <- list(c("ANT1 ME", "IAI ME"), 
                         c("ANT1 ME", "B6 ME"), 
                         c("ANT1 ME", "EC77 ME"),
                         c("IAI ME", "B6 ME"),
                         c("IAI ME", "EC77 ME"),
                         c("B6 ME", "EC77 ME"))
  p <- ggboxplot(
    forbox, x = "strain", y = "covariate", 
    facet.by = "label2", fill = "label2", ylab = ylab) + 
    stat_compare_means(
      comparisons = my_comparisons, 
      label = "p.signif"
    )
  
  ggsave(filename = fname, plot = p, device = 'png', width = 10, height = 8)
}

# find for each variable if there is an association with strain
for(i in 6:ncol(meta)){
  tmp <- meta
  covar <- colnames(tmp)[i]
  colnames(tmp)[i]  <- 'covariate'
  fit1 <- lm(covariate ~ label2, data = tmp)
  group.rSquared <- format(summary(fit1)$r.squared, digits = 2, scientific = T)
  group.pVal <- format(anova(fit1)$'Pr(>F)'[1], digits = 2, scientific = T)
  if(length(which(!is.na(tmp$covariate))) < 40){
    next
  }
  fit2 <- lm(covariate ~ strain, data = tmp)
  strain.rSquared <- format(summary(fit2)$r.squared, digits = 2, scientific = T)
  strain.pVal <- format(anova(fit2)$'Pr(>F)'[1], digits = 2, scientific = T)
  if(i == 6){
    res <- data.frame(covar, group.rSquared, group.pVal, strain.rSquared, strain.pVal)
  } else {
    tt <- data.frame(covar, group.rSquared, group.pVal, strain.rSquared, strain.pVal)
    res <- rbind(res, tt)
  }
}
write.table(res, file = 'results/covars/covar_association.txt', quote = F, row.names = F, sep = "\t")

# correlation with gene exp
# test <- meta %>%
#   mutate(label = label2) %>%
#   dplyr::select(sample, label, strain, weight, VO2max, running.time) %>%
#   cbind(gene = tmp$gene)
# # dput(test)
# summary(lm(gene ~ VO2max + label, data = test))
# summary(lm(gene ~ VO2max + label + strain, data = test))
# summary(lm(gene ~ label + strain, data = test))
# fit <- lm(gene ~ VO2max + label +  strain, data = test)
# res <- tidy(fit)
# res %>%
#   top
# rSquared <- summary(fit)$r.squared
# pVal <- anova(fit)$'Pr(>F)'[1]

