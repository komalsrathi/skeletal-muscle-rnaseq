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
# apply function on each strain
diff.expr <- function(expr, meta, type = "non_exercised", var = "label", batch.correct = NULL, fc = 0, fname = "plot_name", plx = 7){

  strain <- unique(meta$strain)
  print(paste0("Strain:", strain))

  # filter meta file by type
  if(type == "non_exercised"){
    meta <- meta %>%
      mutate(label = ifelse(label == "non_exercised", "non_exercised", "exercised"))
  } else {
    meta <- meta %>%
      filter(label != "non_exercised")
  }


  # filter gene expression matrix by samples
  rownames(meta) <- meta$sample
  expr <- expr[,rownames(meta)]

  # filter gene expression for low expression
  keep.exprs <- filterByExpr(expr)
  expr <- expr[keep.exprs,]

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
  # print(dim(design))

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
  pca.fname <- paste0(strain,'-',fname, '.pdf')
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
    print(comp)
    comp <- gsub(' ','_',comp)
    outputLimma <- topTable(fit2, coef = i, number = Inf)

    # plot volcano
    volc.fname <- paste0(strain,'-',fname, '-volcano.pdf')
    if(!dir.exists('results/volcano')){
      system('mkdir -p results/volcano')
    }
    plotVolcano(result = outputLimma, fname = file.path('results/volcano', volc.fname), yaxis = "P.Value", lfcutoff = 0, pvalcutoff = 0.05)

    if(nrow(outputLimma) > 0){
      outputLimma[,paste0(comp,'_logFC')] <- outputLimma$logFC
      rev.comp <- paste(rev(strsplit(comp, "-")[[1]]), collapse = "-")
      outputLimma[,paste0(rev.comp,'_logFC')] <- -1*(outputLimma$logFC)

      # annotate and filter output
      outputLimma <- annotate.limma(x = outputLimma, foldchange = fc)
      outputLimma$strain <- strain
      print(paste0("Output dimensions:", dim(outputLimma)))
    } else {
      print("No significant difference")
    }
    # newList[[i]] <- assign(comp, outputLimma)
  }

  # return object
  return(assign(comp, outputLimma))
}

# meta  <- meta %>%
#   filter(strain == "ANT1 ME")
# type = "non-exercised"
# expr <- expr.counts[[1]]
# var = "label"
# apply function on each strain
# within strain comparison: Exercised combined (responders + non-responders) vs Non-exercised
non.exercised.vs.exercised <- plyr::dlply(.data = meta, .variables = 'strain',
            .fun = function(x)  diff.expr(expr = expr.counts.mat,
                                          meta = x,
                                          var = "label",
                                          type = "non_exercised",
                                          fname = "within-strain-ex-vs-nonex", plx = 3))

# write to excel file
system('mkdir -p results/summary')
xls.fname <- file.path('results/summary/within-strain-ex-vs-nonex.xlsx')
for (i in 1:length(non.exercised.vs.exercised)){
  if(nrow(non.exercised.vs.exercised[[i]]) > 1){
    sheetname <- unique(as.character(non.exercised.vs.exercised[[i]]$strain))
    write.xlsx(x = non.exercised.vs.exercised[[i]], file = xls.fname, row.names = FALSE, sheetName = sheetname, append=T)
  }
}

# within strain comparison: Exercised (responders) vs Exercised (non-responders)
exercisedres.vs.exercisednores <- plyr::dlply(.data = meta, .variables = 'strain',
                   .fun = function(x)  diff.expr(expr = expr.counts.mat,
                                                 meta = x,
                                                 var = "label",
                                                 type = "exercised",
                                                 fname = "within-strain-exres-vs-exnonres", plx = 1))
# write to excel file
system('mkdir -p results/summary')
xls.fname <- file.path('results/summary/within-strain-exres-vs-exnonres.xlsx')
for (i in 1:length(exercisedres.vs.exercisednores)){
  if(nrow(exercisedres.vs.exercisednores[[i]]) > 1){
    sheetname <- unique(as.character(exercisedres.vs.exercisednores[[i]]$strain))
    write.xlsx(x = exercisedres.vs.exercisednores[[i]], file = xls.fname, row.names = FALSE, sheetName = sheetname, append=T)
  }
}
