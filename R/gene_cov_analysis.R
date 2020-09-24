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
              help = "Metadata file for samples (.txt)"),
  make_option(c("--covariate_file"), type = "character",
              help = "file with covariate data (.txt)"),
  make_option(c("--var_filter"), type = "character", default = TRUE,
              help = "Variance filter: T or F"),
  make_option(c("--col"), type = "character",
              help = "column for comparison"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
counts_matrix <- opt$counts_matrix
meta_file <- opt$meta_file
covariate_file <- opt$covariate_file
var_filter <- opt$var_filter
col <- opt$col
outdir <- opt$outdir

# read data
load(counts_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim(meta_file, stringsAsFactors = F)
rownames(meta) <- meta$sample
meta$label <- ifelse(meta$label == "non_exercised", meta$label, "exercised")
expr.counts.mat <- expr.counts.mat[,rownames(meta)]

# covariates
covariates <- read.delim(covariate_file)
covariates <-  covariates %>%
  column_to_rownames("mice")
covariates <- covariates[rownames(meta),]

# merge covariates to metadata
if(identical(rownames(covariates), rownames(meta))){
  print("identical")
  meta <- cbind(meta, covariates)
}

# function
reg.analysis <- function(x, df){
  x <- x %>%
    inner_join(df %>% rownames_to_column("sample"), by = 'sample')
  covar <- colnames(x)[ncol(x)]
  colnames(x)[ncol(x)] <- 'covar'
  
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

covar_analysis <- function(expr, var = col, var_filter = TRUE, covariates, outdir){
  
  # create output directories
  dir.create(file.path(outdir, 'plots'), showWarnings = F, recursive = T)
  
  # filter expression
  print(var_filter)
  expr <- filterExpr(expr, var.filter = var_filter)
  
  # create design
  print('create design')
  var <- factor(meta[,var])
  design <- model.matrix(~0+var)
  colnames(design) <- levels(var)
  rownames(design) <- meta$sample
  print(dim(design))
  
  # voom normalize data
  print('voom normalization')
  y <- DGEList(counts = as.matrix(expr), genes = rownames(expr))
  y <- calcNormFactors(y)
  v <- voom(counts = y, design = design, plot = FALSE)
  voomData <- v$E
  
  # transpose matrix
  print('transpose')
  voomData.trn <- melt(as.matrix(voomData), varnames = c("gene","sample"), value.name = "expr")
  
  # correlation with gene exp (overall)
  print("Step 1...")
  batch.col <- grep('batch', colnames(meta))
  n <- batch.col + 1
  print(n)
  covars <- colnames(meta[,n:ncol(meta)])
  for(i in 1:length(covars)){
    covar <- covars[i]
    fname <- file.path(outdir, paste0(covar, '_model_output.xlsx'))
    if(!file.exists(fname)){
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
      write.xlsx(x = res, file = fname, row.names = FALSE)
    } else {
      print(fname)
      print("file exists")
    }
  }
  
  # repeat this for each strain
  print("Step 2...")
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
    
    covars <- colnames(subset.meta[,n:ncol(subset.meta)])
    for(i in 1:length(covars)){
      covar <- covars[i]
      print(covar)
      fname <- file.path(outdir, paste0(st, '_', covar, '_model_output.xlsx'))
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
      } else {
        print(fname)
        print("file exists")
      }
    }
  }
  
  # find for each variable if there is an association with strain
  print("Step 3...")
  for(i in n:ncol(meta)){
    tmp <- meta
    covar <- colnames(tmp)[i]
    colnames(tmp)[i]  <- 'covariate'
    fit1 <- lm(covariate ~ label, data = tmp)
    group.rSquared <- format(summary(fit1)$r.squared, digits = 2, scientific = T)
    group.pVal <- format(anova(fit1)$'Pr(>F)'[1], digits = 2, scientific = T)
    if(length(which(!is.na(tmp$covariate))) < 40){
      next
    }
    fit2 <- lm(covariate ~ strain, data = tmp)
    strain.rSquared <- format(summary(fit2)$r.squared, digits = 2, scientific = T)
    strain.pVal <- format(anova(fit2)$'Pr(>F)'[1], digits = 2, scientific = T)
    if(i == n){
      res <- data.frame(covar, group.rSquared, group.pVal, strain.rSquared, strain.pVal)
    } else {
      tt <- data.frame(covar, group.rSquared, group.pVal, strain.rSquared, strain.pVal)
      res <- rbind(res, tt)
    }
  }
  write.xlsx(res, file = file.path(outdir, paste0('covar_association_withStrain.xlsx')), row.names = FALSE)
  
  # boxplots of covariates
  print("Step 4...")
  for(i in n:ncol(meta)){
    forbox <- meta
    ylab <- colnames(forbox)[i]
    colnames(forbox)[i]  <- 'covariate'
    fname <- file.path(outdir, 'plots', paste0(ylab, '.png'))
    p <- ggboxplot(
      forbox, x = "strain", y = "covariate", 
      facet.by = "label", fill = "label", ylab = ylab) + 
      stat_compare_means(ref.group = ".all.", label = "p.signif")
    ggsave(filename = fname, plot = p, device = 'png', width = 10, height = 8)
  }
}

# run function
covar_analysis(expr = expr.counts.mat, 
               var = col, var_filter = var_filter, 
               covariates = covariates, outdir = outdir)






