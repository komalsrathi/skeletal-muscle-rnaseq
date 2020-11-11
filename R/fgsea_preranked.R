# Author: Komal S. Rathi
# Function: Script to generate input lists for GSEA preranked
setwd('~/Projects/Skeletal_Muscle_Patrick/')

# using GSEA preranked: 
# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results

library(readxl)
library(writexl)
library(tidyverse)
library(biomaRt)
library(optparse)
library(msigdbr)
library(fgsea)

option_list <- list(
  make_option(c("--input_dir"), type = "character",
              help = "Input directory with all files"),
  make_option(c("--covars"), type = "character",
              help = "comma separated covariables"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory to write output")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
input_dir <- opt$input_dir
output_dir <- opt$output_dir
covars <- opt$covars
covars <- trimws(strsplit(covars,",")[[1]]) 
covars <- file.path(input_dir, paste0(covars, "_model_output.xlsx"))

# create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# function for fgsea preranked analysis
fgsea_preranked <- function(covar_filename) {
  
  # read covariable output
  dat <- readxl::read_excel(path = covar_filename, sheet = 1)
  res <- dat %>%
    dplyr::select(gene, estimate, model.rSquared, p.value)
  
  # ranked list
  pos <- res %>%
    filter(estimate > 0) %>%
    arrange(p.value, desc(model.rSquared))
  pos <- pos %>%
    mutate(rank = seq(from = 1, to = nrow(pos), by = 1))
  neg <- res %>%
    filter(estimate < 0) %>%
    arrange(p.value, desc(model.rSquared)) 
  neg <- neg %>%
    mutate(rank = rev(seq(from = (nrow(pos) + 1), to = nrow(res), by = 1)))
  res <- rbind(pos, neg)
  res <- res %>%
    arrange(rank)
  fname <- gsub('.*/|_model_output.xlsx', '', covar_filename)
  fname <- file.path(output_dir, paste0(fname, '_ranked_list.xlsx'))
  
  # save ranked list
  print(fname)
  writexl::write_xlsx(x = res, path = fname, format_headers = TRUE)

  # perform fgsea
  # ranks
  ranks <- res$rank
  names(ranks) <- res$gene
  
  # pathways
  h_df = msigdbr(species = "Mus musculus", category = c("H"))
  h_list = h_df %>% split(x = .$gene_symbol, f = .$gs_name)
  c2_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
  c2_list = c2_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pathways = c(h_list, c2_list)
  
  # run fgsea
  fgseaRes <- fgsea(pathways = pathways, 
                    stats    = ranks,
                    eps      = 0.0,
                    minSize  = 15,
                    maxSize  = 1500)
  fgseaRes <- fgseaRes %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(genes = sapply(leadingEdge, paste, collapse=","),
           direction = ifelse(ES > 0, "up", "down"))

  # save output
  fname <- gsub('ranked_list', 'fgsea_output', fname)
  print(fname)
  writexl::write_xlsx(x = fgseaRes, path = fname, format_headers = TRUE)
}

lapply(X = covars, FUN = function(x) fgsea_preranked(covar_filename = x))
