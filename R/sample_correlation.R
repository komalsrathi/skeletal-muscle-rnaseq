# Author: Komal S. Rathi
# Date: 04/20/2020
# Function: Sample correlation heatmap

library(pheatmap)
library(edgeR)

load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$litter[meta$litter == ""] <- "U"
rownames(meta) <- meta$sample
meta$label2 <- ifelse(meta$label == "non_exercised", meta$label, "exercised")
expr.counts.mat <- expr.counts.mat[,rownames(meta)]

# filter low expression
keep.exprs <- filterByExpr(expr.counts.mat)
expr.counts.mat <- expr.counts.mat[keep.exprs,]

# calculate correlation
myCor <- cor(log2(expr.counts.mat+1), method = "pearson")
myCorNames <- rownames(myCor)
pdf(file = 'results/qc/sample_correlation.pdf', width = 20, height = 15)
pheatmap(myCor, annotation_row = meta[c("label2","label","strain")],
         annotation_col = meta[c("label2","label","strain")],
         main = "Sample Correlation Heatmap\n",fontsize = 12,
         display_numbers = TRUE)
dev.off()

