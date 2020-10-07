# Author: Komal S. Rathi
# Function: Batch correction using sva::ComBat

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sva))

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# parameters
option_list <- list(
  make_option(c("--combined_mat"), type = "character",
              help = "Combined expression matrix with multiple batches (Counts) (RData)"),
  make_option(c("--combined_clin"), type = "character",
              help = "Combined clinical file with multiple batches (.tsv)"),
  make_option(c("--corrected_outfile"), type = "character",
              help = "Output filename (RData)"),
  make_option(c("--type"), type = "character", default = 'counts',
              help = "Type of input (counts or fpkm)"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path"))

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
combined_mat <- opt$combined_mat
combined_clin <- opt$combined_clin
corrected_outfile <- opt$corrected_outfile
outdir <- opt$outdir
corrected_outfile <- file.path(outdir, corrected_outfile)
type <- opt$type

# read input data
combined_mat <- get(load(combined_mat))
annot <- combined_mat[[2]]
combined_mat <- combined_mat[[1]]
combined_clin <- read.delim(combined_clin, stringsAsFactors = F)
combined_clin <- combined_clin %>%
  mutate(tmp = sample) %>%
  column_to_rownames('tmp')

# arrange sample identifiers using metadata
combined_mat <- combined_mat[,rownames(combined_clin)]

if(identical(rownames(combined_clin), colnames(combined_mat))){
  print("Matching dimensions")
} else {
  print("Check inputs")
  break
}

# batch correct using ComBat (log2(counts + 1))
# corrected_mat <- ComBat(dat = log2(combined_mat + 1), batch = combined_clin$batch)
# corrected_mat <- 2^(corrected_mat) # back-transform

# batch correct using combat_seq
batch <- combined_clin$batch
cov1 <- as.numeric(factor(combined_clin$treat))
cov2 <- as.numeric(factor(combined_clin$label))
covar_mat <- cbind(cov1, cov2)
corrected_mat <- ComBat_seq(counts = as.matrix(combined_mat), 
                           batch = batch, 
                           group = NULL, 
                           covar_mod = covar_mat)

# save corrected mat and annotation
if(type == "counts"){
  expr.counts <- list(corrected_mat, annot)
  save(expr.counts, file = corrected_outfile)
}  else {
  expr.fpkm <- list(corrected_mat, annot)
  save(expr.fpkm, file = corrected_outfile)
}
