# Author: Komal S. Rathi
# Function: Gene expression correlation with covariates analysis
# For later

library(usedist)
library(edgeR)
library(tidyverse)
library(xlsx)
library(reshape2)
library(bioplotr)
library(broom)
library(ggpubr)
library(optparse)

source('R/filterExpr.R')
source('../Utils/pubTheme.R')

option_list <- list(
  make_option(c("--counts_matrix"), type = "character",
              help = "RData object of counts"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file for samples (.tsv)"),
  make_option(c("--covariate_file"), type = "character",
              help = "file with covariate data"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
counts_matrix <- opt$counts_matrix
meta_file <- opt$meta_file

# read data
load(counts_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim(meta_file, stringsAsFactors = F)
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
voomData.trn <- melt(as.matrix(voomData), varnames = c("gene","sample"), value.name = "expr")

# function
reg.analysis <- function(x, df){
  x <- x %>%
    inner_join(df %>% rownames_to_column("sample"), by = 'sample')
  covar <- colnames(x)[4]
  colnames(x)[4] <- 'covar'
  
  # fit model
  fit <- lm(expr ~ covar, data = x)
  
  # extract info 
  res <- tidy(fit)
  res <- res %>%
    filter(term != "(Intercept)") %>%
    mutate(term = covar) %>%
    as.data.frame()
  
  # add overall model p-value and r-squared
  res$model.rSquared <- summary(fit)$r.squared
  # res$model.pValue <- anova(fit)$'Pr(>F)'[1] same as covariate pvalue
  
  return(res)
}

# correlation with gene exp (overall)
covars <- colnames(meta[,6:ncol(meta)])
for(i in 1:length(covars)){
  covar <- covars[i]
  print(covar)
  df <- meta %>%
    dplyr::select(covar)
  res <- plyr::ddply(voomData.trn, .variables = 'gene', .fun = function(x) reg.analysis(x, df))
  
  # adjust pvalues for each
  res <- res %>%
    mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    mutate(correlated = ifelse(fdr < 0.05, TRUE, FALSE)) %>%
    arrange(fdr) 
  
  # write out files
  fname <- paste0('results/covars/',covar,'_model_output.xlsx')
  write.xlsx(x = res, file = fname, row.names = FALSE)
}

# repeat this for each strain
strains <- unique(meta$strain)
for(j in 1:length(strains)){
  print(strains[j])
  st <- gsub(' ','_',strains[j])
  subset.meta <- meta %>%
    mutate(tmp = sample) %>%
    filter(strain == strains[j]) %>%
    column_to_rownames('tmp')
  
  subset.trn <- voomData.trn %>%
    filter(sample %in% subset.meta$sample)
  
  covars <- colnames(subset.meta[,6:ncol(subset.meta)])
  for(i in 1:length(covars)){
    covar <- covars[i]
    print(covar)
    fname <- paste0('results/covars/',st,'_',covar,'_model_output.xlsx')
    if(!file.exists(fname)){
      df <- subset.meta %>%
        dplyr::select(covar)
      if(all(is.na(df[,covar]))){
        next
      }
      res <- plyr::ddply(subset.trn, .variables = 'gene', .fun = function(x) reg.analysis(x, df))
      
      # adjust pvalues for each
      res <- res %>%
        mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
        mutate(correlated = ifelse(fdr < 0.05, TRUE, FALSE)) %>%
        arrange(fdr) 
      
      # write out files
      write.xlsx(x = res, file = fname, row.names = FALSE)
    }
  }
}

# boxplots of covariates
for(i in 6:ncol(meta)){
  forbox <- meta
  ylab <- colnames(forbox)[i]
  colnames(forbox)[i]  <- 'covariate'
  fname <- paste0('results/covars/plots/',ylab,'.png')
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
write.xlsx(res, file = 'results/covars/covar_association_withStrain.xlsx', row.names = FALSE)
