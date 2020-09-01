# Author: Komal S. Rathi
# Date: 04/20/2020
# Function: Sample correlation heatmap

library(optparse)
library(pheatmap)
library(edgeR)

option_list <- list(
  make_option(c("--counts_matrix"), type = "character",
              help = "RData object of counts"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file for samples (.tsv)"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
counts_matrix <- opt$counts_matrix
meta_file <- opt$meta_file
outdir <- opt$outdir

# make output directories
outdir <- file.path(outdir, 'clustering-qc')
dir.create(outdir, showWarnings = F, recursive = TRUE)

# load data
load(counts_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
meta <- read.delim(meta_file, stringsAsFactors = F)
rownames(meta) <- meta$sample
meta$label <- ifelse(meta$label == "non_exercised", meta$label, "exercised")
expr.counts.mat <- expr.counts.mat[,rownames(meta)]

# filter low expression
keep.exprs <- filterByExpr(expr.counts.mat)
expr.counts.mat <- expr.counts.mat[keep.exprs,]

# calculate correlation
myCor <- cor(log2(expr.counts.mat+1), method = "pearson")
myCorNames <- rownames(myCor)
fname <- file.path(outdir, 'sample_correlation.pdf')
pdf(file = fname, width = 20, height = 15)
pheatmap(myCor, annotation_row = meta[c("label","strain")],
         annotation_col = meta[c("label","strain")],
         main = "Sample Correlation Heatmap\n", fontsize = 12,
         display_numbers = TRUE)
dev.off()

