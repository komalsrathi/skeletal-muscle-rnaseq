# Author: Komal S. Rathi
# Function: To create collapsed matrices of gene symbols and samples
# Date: 04/02/2020

setwd('~/Projects/Skeletal_Muscle_Patrick/')

library(tidyverse)

# load RSEM data
load('data/skeletal_muscle_mm10.RData')

# annotation
annot <- data.table::fread('data/gencode.vM17.annotation.txt')
annot <- annot %>%
  select(gene_id,  gene_symbol, biotype) %>%
  unique()

# merge annotation with expression
collapse.mat <- function(expr, geneAnnot){
  
  # add annotation
  expr <- geneAnnot %>%
    select(gene_id, gene_symbol)  %>%
    inner_join(expr, by = c('gene_id'))
  
  # collapse to gene symbols
  expr <- expr %>%
    dplyr::mutate(means = rowMeans(.[3:50])) %>%
    arrange(desc(means)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    dplyr::select(-c(means)) %>%
    column_to_rownames(var = "gene_symbol")
  
  # subset annotation
  geneAnnot <- geneAnnot %>%
    filter(gene_id %in% expr$gene_id)
  
  # remove gene_id from matrix
  expr <- expr %>%
    select(-c(gene_id))
  
  newList <- list(expr, geneAnnot)
  return(newList)
}

# create collapsed count matrix
expr.counts <- collapse.mat(expr = expr.counts, geneAnnot = annot)
expr.counts.annot <- expr.counts[[2]]
expr.counts.mat <- expr.counts[[1]]

# create collapsed fpkm matrix
expr.fpkm <- collapse.mat(expr = expr.fpkm, geneAnnot = annot)
expr.fpkm.annot <- expr.fpkm[[2]]
expr.fpkm.mat <- expr.fpkm[[1]]

save(expr.counts, file = "data/collapsed_counts_matrix.RData")
save(expr.fpkm, file = "data/collapsed_fpkm_matrix.RData")
