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

# Analysis:
# 1. Between strain comparison: Non-exercised mice
# 2. Between strain comparison: Exercised mice combined (responders + non-responders)

# generalized function
diff.expr <- function(expr, meta, annot, type = "non_exercised", var = "strain", batch.correct = NULL, fc = 1, fname = "plot_name", plx = 7, write_to_excel = FALSE, write_to_text = TRUE){

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

  # batch correction
  if(!is.null(batch.correct)){
    batch <- factor(meta[,batch.correct])
    voomData <- ComBat(voomData, batch = batch, mod = design, par.prior = T)
  }

  # tsne and pca plot of voom normalized data
  if(!dir.exists('results/pca')){
    system('mkdir -p results/pca')
  }
  pca.fname <- paste0(fname, '.pdf')
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
  newList <- list()
  for(i in 1:t){
    comp <- colnames(fit2$coefficients)[i] # comparison
    comp <- gsub(' ','_',comp)
    print(comp)
    outputLimma <- topTable(fit2, coef = i, number = Inf)

    # plot volcano
    volc.fname <- paste0(comp,'-',fname, '-volcano.pdf')
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
      outputLimma$comparison <- comp
      print(paste0("Output dimensions:", dim(outputLimma)))
    } else {
      outputLimma <- data.frame(gene_symbol = NULL)
      print("No significant difference")
    }
    newList[[i]] <- assign(comp, outputLimma)
  }

  # write to excel file
  if(write_to_excel == TRUE){
    system('mkdir -p results/summary')
    xls.fname <- file.path('results/summary',paste0(fname,'.xlsx'))
    for (i in 1:length(newList)){
      if(nrow(newList[[i]]) > 1){
        sheetname <- unique(as.character(newList[[i]]$comp))
        write.xlsx(x = newList[i], file = xls.fname, row.names = FALSE, sheetName = sheetname, append=T)
        gc()
      }
    }
  }
  
  # write to text file
  if(write_to_text == TRUE){
    system('mkdir -p results/summary')
    for(i in 1:length(newList)){
      if(nrow(newList[[i]] > 1)){
        sheetname <- unique(as.character(newList[[i]]$comp))
        txt.fname <- file.path('results/summary', paste0(fname, '-',sheetname,'.txt'))
        tmp <- newList[[i]] %>%
          filter(adj.P.Val < 0.05)
        if(nrow(tmp) > 0) {
          write.table(tmp, file = txt.fname, quote = F, row.names = F, sep = "\t")
        }
      }
    }
  }

  # return object
  return(newList)
}

# 1. between strain comparison of non-exercised mice
nonex.between.strains <- diff.expr(expr = expr.counts[[1]], meta = meta,
                                   annot = expr.counts.annot,
                                   type = "non_exercised",
                                   var = "strain", batch.correct = NULL, fc = 0, 
                                   fname = "between-strain-nonex", plx = 7, 
                                   write_to_excel = FALSE, 
                                   write_to_text = TRUE)

# 2. between strain comparison of exercised mice (responders + non-responders)
ex.between.strains <- diff.expr(expr = expr.counts[[1]], meta = meta,
                                annot = expr.counts.annot,
                                type = c("exercised_responders","exercised_non_responders"),
                                var = "strain", batch.correct = NULL, 
                                fc = 0, fname = "between-strain-ex", plx = 7, 
                                write_to_excel = FALSE,
                                write_to_text = TRUE)

# qc: check differences between exercised responders
ex.res.between.strains <- diff.expr(expr = expr.counts[[1]], meta = meta,
                                    annot = expr.counts.annot,
                                    type = c("exercised_responders"),
                                    var = "strain", batch.correct = NULL, 
                                    fc = 0, fname = "between-strain-exres", plx = 3, 
                                    write_to_excel = FALSE, 
                                    write_to_text = TRUE)

# qc: check differences between exercised non-responders
ex.nonres.between.strains <- diff.expr(expr = expr.counts[[1]], meta = meta,
                                       annot = expr.counts.annot,
                                       type = c("exercised_non_responders"),
                                       var = "strain", batch.correct = NULL, fc = 0, 
                                       fname = "between-strain-exnonres", plx = 3, 
                                       write_to_excel = FALSE,
                                       write_to_text = TRUE)

# observation:
# in all four cases, differences are seen between A vs E, I and B but not between other strains.
