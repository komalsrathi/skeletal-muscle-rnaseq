# This script will filter genes using a bunch of metrics
# 1. Low Expression
# 2. IQR: Low Variance
library(genefilter)
library(edgeR)

filterExpr <- function(expr.counts.mat, var.filter) {
  
  # 1. filter by expression
  print('Filtering by expression...')
  keep.exprs <- filterByExpr(expr.counts.mat, min.count = 10)
  expr.counts.mat <- expr.counts.mat[keep.exprs,]
  print(dim(expr.counts.mat))
  
  # 2. filter by IQR (low variance genes)
  if(var.filter){
    print('Filtering by variance...')
    expr.counts.mat <- varFilter(as.matrix(expr.counts.mat), 
                                 var.func = IQR, 
                                 var.cutoff = 0.5, 
                                 filterByQuantile = TRUE)
    print(dim(expr.counts.mat))
  }
  
  # return filtered matrix
  return(expr.counts.mat)
}


