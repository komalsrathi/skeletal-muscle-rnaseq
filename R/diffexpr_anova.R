# Author: Komal S. Rathi
# Date: 04/30/2020
# Function: ANOVA across 4 strains * 2 types
library(reshape2)
library(pheatmap)
library(xlsx)

# read data
load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$strain <- gsub(' ','_',meta$strain)
meta$litter[meta$litter == ""] <- "U"
meta$label <- ifelse(meta$label == "non_exercised", meta$label, "exercised")

# define 8 groups
meta$group <- paste0(meta$strain, '-',meta$label)
# expr <- expr.counts.mat
# lb <- 'exercised'

calc.var <- function(expr, meta, lb = NULL, fname){

  # filter by label
  if(!is.null(lb)) {
    meta <- meta %>%
      filter(label == lb)
  }

  # filter gene expression for low expression
  rownames(meta) <- meta$sample
  expr <- expr[,rownames(meta)]
  keep.exprs <- filterByExpr(expr)
  expr <- expr[keep.exprs,]

  if(identical(rownames(meta), colnames(expr))) {
    print("Proceed")
  } else {
    break
  }
  expr <- melt(as.matrix(expr))

  # for each gene in the matrix, do anova
  aov.test <- function(myDataExp, myMeta) {
    myMeta <- myMeta %>%
      inner_join(myDataExp, by = c("sample" = "Var2"))
    anovaRes <- anova(lm(value ~ group, data = myMeta))
    res <- anovaRes[1,5]
    return(res)
  }
  aov.res <- plyr::ddply(.data = expr, .variables = 'Var1', .fun = function(x)  aov.test(myDataExp = x, myMeta = meta))
  colnames(aov.res) <- c("Gene","Pval")

  # adjust p-value using Bonferroni
  aov.res$AdjPval <- p.adjust(aov.res$Pval, method = "bonf")
  aov.res <- aov.res %>%
    mutate(DEGAnnot = ifelse(AdjPval < 0.05, TRUE, FALSE))
  if(!dir.exists('results/anova')){
    system('mkdir -p results/anova')
  }
  table.fname <- file.path('results/anova/', paste0(fname, '.xlsx'))
  write.xlsx(x = aov.res, file = table.fname, row.names = FALSE, sheetName = lb, append=T)

  # heatmap of adj pvalue < 0.05
  for.heatmap <- aov.res %>%
    filter(DEGAnnot) %>%
    .$Gene %>%
    as.character()

  # FPKM to z-score
  load('data/collapsed_fpkm_matrix.RData')
  fpkm.mat <- expr.fpkm[[1]]
  fpkm.mat <- fpkm.mat[for.heatmap,rownames(meta)]

  # function to calculate z-score
  getZ <- function(x) {
    x <- log2(x+1)
    out <- (x-mean(x))/sd(x)
    return(out)
  }
  fpkm.mat <- as.data.frame(t(apply(fpkm.mat, FUN = getZ, MARGIN = 1)))
  fpkm.meta <- meta[colnames(fpkm.mat),c("strain","label")]

  # plot pheatmap
  heatmap.fname <- file.path('results/anova/', paste0(fname, '_heatmap.pdf'))
  pdf(file = heatmap.fname, width = 10, height = 15)
  pheatmap(fpkm.mat, cluster_cols = FALSE,
           main = paste("Genes (ANOVA P-adj < 0.05)", ":", nrow(fpkm.mat)),
           fontsize_col = 6, fontsize_row = 6,
           annotation_col = fpkm.meta, scale = "row", angle_col = 45)
  dev.off()
}

calc.var(expr = expr.counts.mat, meta = meta, lb = 'exercised', fname = 'exercised-anova')
calc.var(expr = expr.counts.mat, meta = meta, lb = 'non_exercised', fname = 'non-exercised-anova')

