# Author: Komal S. Rathi
# Differential expression (Human)
# Date: 09/02/2020

options(java.parameters = "- Xmx1024m")
library(optparse)
library(tidyverse)
library(limma)
library(edgeR)
library(sva)
library(xlsx)

option_list <- list(
  make_option(c("--counts_matrix"), type = "character",
              help = "RData object of counts"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file for samples (.tsv)"),
  make_option(c("--type"), type = "character",
              help = "Type of comparison"),
  make_option(c("--col"), type = "character",
              help = "column for comparison"),
  make_option(c("--fc"), type = "character", default = 0,
              help = "fold-change for output files"),
  make_option(c("--plx"), type = "character", default = 0,
              help = "perplexity for t-SNE"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files"),
  make_option(c("--excel"), type = "character",
              default = TRUE,
              help = "Write to excel file?"),
  make_option(c("--text"), type = "character",
              default = TRUE,
              help = "Write to txt file?"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
counts_matrix <- opt$counts_matrix
meta_file <- opt$meta_file
type <- opt$type
col <- opt$col
fc <- opt$fc
plx <- opt$plx
prefix <- opt$prefix
excel <- opt$excel
text <- opt$text
outdir <- opt$outdir

# source functions
source('~/Projects/Utils/design_pairs.R')
source('R/annotate_limma.R')
source('R/PCA-plot.R')
source('R/Volcano-plot.R')
source('R/filterExpr.R')

# read data
load(counts_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim(meta_file, stringsAsFactors = F)
type <- trimws(strsplit(type,",")[[1]]) 
plx <- as.numeric(plx)
fc <- as.numeric(fc)

# make output directories
pca.out <- file.path(outdir, 'pca')
dir.create(pca.out, showWarnings = F, recursive = TRUE)
tsne.out <- file.path(outdir, 'tsne')
dir.create(tsne.out, showWarnings = F, recursive = TRUE)
volcano.out <- file.path(outdir, 'volcano')
dir.create(volcano.out, showWarnings = F, recursive = TRUE)
summary.out <- file.path(outdir, 'summary')
dir.create(summary.out, showWarnings = F, recursive = TRUE)

# Analysis:
# 1. All comparisons between cell lines

# generalized function
diff.expr <- function(expr, meta, annot, type = c("exercised_responders","exercised_non_responders"), 
                      var = 'label', fc = 0, plx = 7, fname = "plot_name", 
                      write_to_excel = TRUE,  write_to_text = TRUE){
  
  # filter meta file by type
  meta <- meta %>%
    filter(label %in% type) %>%
    arrange(sample)
  
  # filter gene expression matrix by samples
  rownames(meta) <- meta$sample
  expr <- expr[,rownames(meta)]
  
  # filter gene expression for low expression
  # keep.exprs <- filterByExpr(expr)
  # expr <- expr[keep.exprs,]
  expr <- filterExpr(expr.counts.mat = expr)
  
  if(identical(rownames(meta), colnames(expr))) {
    print("Proceed")
  } else {
    break
  }
  
  # create design
  var <- factor(meta[,var])
  design <- model.matrix(~0+var)
  colnames(design) <- levels(var)
  rownames(design) <- meta$sample
  print(dim(design))
  
  # voom normalize
  y <- DGEList(counts = as.matrix(expr), genes = rownames(expr))
  y <- calcNormFactors(y)
  v <- voom(counts = y, design = design, plot = FALSE)
  voomData <- v$E
  
  # tsne and pca plot of voom normalized data
  pca.fname <- paste0(fname, '.pdf')
  pca.plot(voomData = voomData, meta = meta, fname = file.path(pca.out, pca.fname), color_var = 'label', shape_var = 'label')
  tsne.plot(voomData = voomData, meta = meta, fname = file.path(tsne.out, pca.fname), plx = plx, color_var = 'label', shape_var = 'label')
  
  # all levels of design
  all.pairs <- design.pairs(levels = colnames(design))
  
  # differential expression
  fit <- lmFit(voomData, design)
  fit2 <- contrasts.fit(fit, all.pairs)
  fit2 <- eBayes(fit2)
  
  # create list of all combinations
  t <- ncol(fit2$coefficients)
  newList <- list()
  for(i in 1:t){
    comp <- colnames(fit2$coefficients)[i] # comparison
    comp <- gsub(' ','_',comp)
    outputLimma <- topTable(fit2, coef = i, number = Inf)
    
    # plot volcano
    volc.fname <- paste0(fname, '-volcano.pdf')
    plotVolcano(result = outputLimma, fname = file.path(volcano.out, volc.fname), yaxis = "P.Value", lfcutoff = 0, pvalcutoff = 0.05)
    
    if(nrow(outputLimma) > 0){
      outputLimma[,paste0(comp,'_logFC')] <- outputLimma$logFC
      rev.comp <- paste(rev(strsplit(comp, "-")[[1]]), collapse = "-")
      outputLimma[,paste0(rev.comp,'_logFC')] <- -1*(outputLimma$logFC)
      
      # annotate and filter output
      outputLimma <- annotate.limma(x = outputLimma, foldchange = fc, annot)
      print(paste0("Output dimensions:", dim(outputLimma)))
    } else {
      print("No significant difference")
    }
    newList[[i]] <- outputLimma
  }
  
  # write to excel file
  if(write_to_excel == TRUE){
    print("Writing output to excel...")
    xls.fname <- file.path(summary.out, paste0(fname,'.xlsx'))
    for (i in 1:length(newList)){
      if(nrow(newList[[i]]) > 1){
        write.xlsx(x = newList[[i]], file = xls.fname, row.names = FALSE, sheetName = "all_comparisons", append=T)
        gc()
      }
    }
  }
  
  # write to text file
  if(write_to_text == TRUE){
    print("Writing output to text...")
    txt.fname <- file.path(summary.out, paste0(fname, '.txt'))
    for(i in 1:length(newList)){
      if(nrow(newList[[i]]) > 1){
        tmp <- newList[[i]] %>%
          filter(adj.P.Val < 0.05)
        if(nrow(tmp) > 0) {
          write.table(tmp, file = txt.fname, quote = F, row.names = F, sep = "\t")
        }
      }
    }
  }
}

# call function
diff.expr(expr = expr.counts.mat, meta = meta, annot = expr.counts.annot, 
          type = type, var = col, fc = fc, plx = plx, 
          fname = prefix, write_to_excel = excel, write_to_text = text)
