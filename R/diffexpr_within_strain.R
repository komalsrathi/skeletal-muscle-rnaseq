# Author: Komal S. Rathi
# Function: Differential expression analysis
# Date: 04/02/2020

options(java.parameters = "- Xmx1024m")
library(tidyverse)
library(limma)
library(edgeR)
library(xlsx)

source('~/Projects/Utils/design_pairs.R')
source('R/annotate_limma.R')
source('R/PCA-plot.R')
source('R/Volcano-plot.R')
source('R/filterExpr.R')

load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$litter[meta$litter == ""] <- "U"

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
  if(!dir.exists('results/pca')){
    system('mkdir -p results/pca')
  }
  pca.fname <- paste0(st,'-',fname, '.pdf')
  pca.plot(voomData = voomData, meta = meta, fname = file.path('results/pca', pca.fname))
  
  if(!dir.exists('results/tsne')){
    system('mkdir -p results/tsne')
  }
  tsne.plot(voomData = voomData, meta = meta, fname = file.path('results/tsne', pca.fname), plx = plx)
  
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
    volc.fname <- paste0(st,'-',fname, '-volcano.pdf')
    if(!dir.exists('results/volcano')){
      system('mkdir -p results/volcano')
    }
    plotVolcano(result = outputLimma, fname = file.path('results/volcano', volc.fname), yaxis = "P.Value", lfcutoff = 0, pvalcutoff = 0.05)
    
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
  system('mkdir -p results/summary')
  
  # write to excel file
  if(write_to_excel == TRUE){
    xls.fname <- file.path('results/summary',paste0(fname,'.xlsx'))
    write.xlsx(x = outputLimma, file = xls.fname, row.names = FALSE, sheetName = strain, append=T)
    gc()
  }
  
  # write to text file
  if(write_to_text == TRUE){
    txt.fname <- file.path('results/summary', paste0(fname,'-',strain, '.txt'))
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
  
  # ex vs nonex
  diff.expr(expr = expr.counts.mat,
            meta = meta, 
            st = st,
            annot = expr.counts.annot,
            var = "label",
            type = "non_exercised", 
            fname = "within-strain-ex-vs-nonex", plx = 3,
            write_to_excel = TRUE, write_to_text = TRUE)
  
  # ex res vs ex non-res
  diff.expr(expr = expr.counts.mat,
            meta = meta,
            st = st,
            annot = expr.counts.annot,
            var = "label",
            type = "exercised",
            fname = "within-strain-exres-vs-exnonres", plx = 1,
            write_to_excel = TRUE, write_to_text = TRUE)
}
