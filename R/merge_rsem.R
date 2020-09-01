# Script to merge RSEM files

library(reshape2)
library(optparse)

option_list <- list(
  make_option(c("--sourcedir"), type = "character",
              help = "Source directory with all RSEM files"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
sourcedir <- opt$sourcedir
prefix <- opt$prefix
outdir <- opt$outdir

# universal function to merge expr and fusion results
merge.res <- function(nm){
  print(nm)
  sample_name <- gsub('.*[/]|[.].*','',nm)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# expression (FPKM)
lf <- list.files(path = sourcedir, pattern = "*.genes.results", recursive = T, full.names = T)
expr.mat <- lapply(lf, FUN = function(x) merge.res(x))
expr.mat <- data.table::rbindlist(expr.mat)
expr.fpkm <- dcast(expr.mat, gene_id~sample_name, value.var = 'FPKM')
expr.counts <- dcast(expr.mat, gene_id~sample_name, value.var = 'expected_count')
output <- file.path(outdir, paste0(prefix, '.RData'))
save(expr.fpkm, expr.counts, file = output)
