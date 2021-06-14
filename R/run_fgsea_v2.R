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
  make_option(c("--type"), type = "character",
              help = "Ascending (a) or Descending (d)"),
  make_option(c("--ranking"), type = "character",
              help = "Ranking metric (p for p-value/r-squared combination and e for estimate)"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory to write output")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
input_dir <- opt$input_dir
output_dir <- opt$output_dir
type <- opt$type
ranking <- opt$ranking
covars <- opt$covars
covars <- trimws(strsplit(covars,",")[[1]]) 

# read files from input directory
lf <- list.files(path = input_dir, full.names = T)
lf <- grep(paste0(covars, collapse = "|"), lf, value = T)
lf <- lf[grep('_exercised', lf)]

# create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# function for fgsea preranked analysis
fgsea_preranked <- function(res, fname){
  
  # perform fgsea
  # ranks
  ranks <- res$rank
  names(ranks) <- res$gene
  
  # pathways
  h_df = msigdbr(species = "Mus musculus", category = c("H"))
  h_list = h_df %>% split(x = .$gene_symbol, f = .$gs_name)
  c2_df = msigdbr(species = "Mus musculus", category = "C2")
  c2_df = c2_df %>%
    dplyr::filter(gs_subcat != "CGP")
  c2_list = c2_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pathways = c(h_list, c2_list)
  
  # run fgsea - classic method
  set.seed(42)
  fgseaRes <- fgsea(pathways = pathways, 
                    gseaParam = 0, 
                    stats = ranks, eps = 0.0,
                    minSize = 15, maxSize = 1500)
  fgseaRes <- fgseaRes %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    mutate(genes = sapply(leadingEdge, paste, collapse=","),
           direction = ifelse(ES > 0, "up", "down"))
  
  # save output
  fname <- gsub('ranked_list', 'fgsea_output', fname)
  fname <- gsub('.rnk', '.xlsx', fname)
  print(fname)
  writexl::write_xlsx(x = fgseaRes, path = fname, format_headers = TRUE)
}

# create ranked list
create_ranked_list <- function(covar_filename) {
  
  # read covariable output
  dat <- readxl::read_excel(path = covar_filename, sheet = 1)
  res <- dat %>%
    dplyr::select(gene, estimate, model.rSquared, p.value)
  
  # ranked list by p-value and r-squared
  if(ranking == "p"){
    print("Ranking by p-value & r-squared")
    var <- "pvalue_rsquared"
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
    # ascending
    res <- res %>%
      arrange(rank)
    # descending
    if(type == "d"){
      print("Descending order")
      res <- res %>% 
        mutate(rank = rev(rank))
    } else {
      print("Ascending order")
    }
  } 
  
  # ranked list by estimate
  if(ranking == "e"){
    print("Ranking by estimate")
    var <- "estimate"
    # res <- res %>%
    #   mutate(rank = estimate)
    if(type == "a"){
      print("Ascending order")
      res <- res %>%
        arrange(estimate) %>%
        mutate(rank = 1:nrow(res))
    } else if(type == "d"){
      print("Descending order")
      res <- res %>%
        arrange(desc(estimate)) %>%
        mutate(rank = 1:nrow(res))
    } 
  }
  
  fname <- gsub('.*/|_model_output_exercised.xlsx', '', covar_filename)
  if(type == "d"){
    fname <- file.path(output_dir, paste(fname, 'ranked_list', var, 'descending_exercised.rnk', sep = "_"))
  } else {
    fname <- file.path(output_dir, paste(fname, 'ranked_list', var, 'ascending_exercised.rnk', sep = "_"))
  }
  
  # save ranked list
  print(fname)
  write.table(x = res, file = fname, sep = "\t", quote = F, col.names = F, row.names = F)
  
  # run fgsea in classic mode
  fgsea_preranked(res = res, fname = fname)
}



lapply(X = lf, FUN = function(x) create_ranked_list(covar_filename = x))
