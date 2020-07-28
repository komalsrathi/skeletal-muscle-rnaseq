# Author: Komal S. Rathi
# Function: Differential expression analysis
# Date: 04/02/2020

options(java.parameters = "- Xmx1024m")
library(tidyverse)
library(limma)
library(edgeR)
library(xlsx)
library(optparse)

option_list <- list(
  make_option(c("--counts_matrix"), type = "character",
              help = "RData object of counts"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file for samples (.tsv)"),
  make_option(c("--type"), type = "character",
              help = "Type of comparison"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
counts_matrix <- opt$counts_matrix
meta_file <- opt$meta_file
type <- opt$type
prefix <- opt$prefix
outdir <- opt$outdir

source('~/Projects/Utils/design_pairs.R')
source('R/annotate_limma.R')
source('R/PCA-plot.R')
source('R/Volcano-plot.R')
source('R/filterExpr.R')

load(counts_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim(meta_file, stringsAsFactors = F)

# Analysis
# 3. Within strain comparison: Exercised combined (responders + non-responders) vs Non-exercised
# 4. Within strain comparison: Exercised (responders) vs Exercised (non-responders)

# generalized function
diff.expr <- function(expr, meta, annot, st, type = "non_exercised", var = "label", batch.correct = NULL, fc = 0, fname = "plot_name", plx = 7, write_to_text = TRUE, write_to_excel = FALSE){
  
  print(paste0("Strain:", st))
  
  # filter meta file by strain and type
  if(type == "non_exercised"){
    meta <- meta %>%
      filter(strain == st) %>%
      mutate(label = ifelse(label == "non_exercised", "non_exercised", "exercised"))
  } else {
    meta <- meta %>%
      filter(strain == st) %>%
      filter(label != "non_exercised")
  }
  
  # filter gene expression matrix by samples
  rownames(meta) <- meta$sample
  expr <- expr[,rownames(meta)]
  
  # filter gene expression for low expression
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
  
  # voom normalize
  y <- DGEList(counts = as.matrix(expr), genes = rownames(expr))
  y <- calcNormFactors(y)
  v <- voom(counts = y, design = design, plot = FALSE)
  voomData <- v$E
  
  # batch correction
  if(!is.null(batch.correct)){
    batch <- factor(exercised[,batch.correct])
    voomData <- ComBat(voomData, batch = batch, mod = design, par.prior = T)
  }
  
  # tsne and pca plot of voom normalized data
  pca.out <- file.path(outdir, 'pca')
  if(!dir.exists(pca.out)){
    system(paste0('mkdir -p ', pca.out))
  }
  pca.fname <- paste0(st,'-',fname, '.pdf')
  pca.plot(voomData = voomData, meta = meta, fname = file.path(pca.out, pca.fname))
  
  tsne.out <- file.path(outdir, 'tsne')
  if(!dir.exists(tsne.out)){
    system(paste0('mkdir -p ', tsne.out))
  }
  tsne.plot(voomData = voomData, meta = meta, fname = file.path(tsne.out, pca.fname), plx = plx)
  
  # all levels of design
  all.pairs <- design.pairs(levels = colnames(design))
  
  # differential expression
  fit <- lmFit(voomData, design)
  fit2 <- contrasts.fit(fit, all.pairs)
  fit2 <- eBayes(fit2)
  
  # create list of all combinations
  t <- ncol(fit2$coefficients)
  for(i in 1:t){
    comp <- colnames(fit2$coefficients)[i] # comparison
    print(comp)
    comp <- gsub(' ','_',comp)
    outputLimma <- topTable(fit2, coef = i, number = Inf)
    
    # plot volcano
    volcano.out <- file.path(outdir, 'volcano')
    volc.fname <- paste0(st,'-',fname, '-volcano.pdf')
    if(!dir.exists(volcano.out)){
      system(paste0('mkdir -p ', volcano.out))
    }
    plotVolcano(result = outputLimma, fname = file.path(volcano.out, volc.fname), yaxis = "P.Value", lfcutoff = 0, pvalcutoff = 0.05)
    
    if(nrow(outputLimma) > 0){
      outputLimma[,paste0(comp,'_logFC')] <- outputLimma$logFC
      rev.comp <- paste(rev(strsplit(comp, "-")[[1]]), collapse = "-")
      outputLimma[,paste0(rev.comp,'_logFC')] <- -1*(outputLimma$logFC)
      
      # annotate and filter output
      outputLimma <- annotate.limma(x = outputLimma, foldchange = fc, annot)
      outputLimma$strain <- st
      print(paste0("Output dimensions:", dim(outputLimma)))
    } else {
      print("No significant difference")
    }
    # newList[[i]] <- assign(comp, outputLimma)
  }
  
  strain  <- gsub(' ','_',st)
  summary.out <- file.path(outdir, 'summary')
  system(paste0('mkdir -p ', summary.out))
  
  # write to excel file
  if(write_to_excel == TRUE){
    xls.fname <- file.path(summary.out, paste0(fname, '.xlsx'))
    write.xlsx(x = outputLimma, file = xls.fname, row.names = FALSE, sheetName = strain, append=T)
    gc()
  }
  
  # write to text file
  if(write_to_text == TRUE){
    txt.fname <- file.path(summary.out, paste0(fname, '-', strain, '.txt'))
    if(nrow(outputLimma) > 1){
      tmp <- outputLimma %>%
        filter(adj.P.Val < 0.05)
      if(nrow(tmp) > 0) {
        write.table(tmp, file = txt.fname, quote = F, row.names = F, sep = "\t")
      }
    }
  }
}

# apply function on each strain
strains <- unique(meta$strain)
for(i in 1:length(strains)){
  st <- strains[i]
  diff.expr(expr = expr.counts.mat,
            meta = meta, 
            st = st,
            annot = expr.counts.annot,
            var = "label",
            type = type, 
            fname = prefix, 
            plx = 3,
            write_to_excel = TRUE, 
            write_to_text = TRUE)
  
}
