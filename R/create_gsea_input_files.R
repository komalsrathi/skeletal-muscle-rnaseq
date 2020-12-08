# Author: Komal S. Rathi
# Function: Create GSEA input files

# libraries
library(tidyverse)
library(reshape2)
library(biomaRt)
library(optparse)

option_list <- list(
  make_option(c("--count_matrix"), type = "character",
              help = "Count matrix (.RData)"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files")
)
# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
count_matrix <- opt$count_matrix
meta_file <- opt$meta_file
output_dir <- opt$output_dir
prefix <- opt$prefix
# output_dir <- 'results/skeletal_muscle/gsea'
dir.create(output_dir, showWarnings = F, recursive = T)

# load count matrix
# count_matrix <- 'data/mouse_skm/skm_collapsed_counts_matrix.RData'
load(count_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.mat <- expr.counts.mat[apply(expr.counts.mat!=0, 1, all),]

# metadata
# meta_file <- 'data/mouse_skm/skm-meta-data.txt'
meta_file <- read.delim(meta_file, stringsAsFactors = F)
meta_file$label <- ifelse(meta_file$label == "non_exercised", "non_exercised", "exercised")
meta_file$strain <- gsub(' ','_',meta_file$strain)

# function to convert mouse to human
convertMouseGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x, 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol","gene_biotype"), 
                   martL = human, uniqueRows=T)
  return(genesV2)
}
human.genes <- convertMouseGeneList(rownames(expr.counts.mat))
human.genes <- human.genes %>%
  filter(HGNC.symbol != "")
expr.counts.mat <- expr.counts.mat %>%
  rownames_to_column("gene_symbol") %>%
  inner_join(human.genes, by = c("gene_symbol" = "MGI.symbol")) %>%
  dplyr::select(-c(gene_symbol, Gene.type)) %>%
  dplyr::mutate(means = rowMeans(dplyr::select(.,-HGNC.symbol))) %>%
  arrange(desc(means)) %>%
  distinct(HGNC.symbol, .keep_all = TRUE) %>% 
  dplyr::select(-c(means)) %>%
  column_to_rownames(var = "HGNC.symbol")

# function
gsea.input <- function(counts_collapsed, meta, groups, strains, type = 'within', gct_file, cls_file) {
  
  # add group to meta file
  meta <- meta %>%
    filter(label %in% groups,
           strain %in% strains) %>%
    mutate(tmp = sample) %>%
    column_to_rownames('tmp') %>%
    arrange(label)
  
  # order matrix
  counts_collapsed <- counts_collapsed[,rownames(meta)]
  
  # print dimensions
  print(strains)
  print(dim(meta))
  print(dim(counts_collapsed))
  
  # gct file
  gct <- counts_collapsed
  add <- data.frame(NAME = c("#1.2", nrow(gct), "NAME"), 
                    Description = c('', ncol(gct), "Description"))
  total.cols <- ncol(gct) + 2
  add[,3:total.cols] <- ''
  colnames(add)[3:total.cols] <- colnames(gct)
  add[3, 3:total.cols] <- colnames(gct)
  annot <- data.frame(NAME = rownames(gct), Description = 'na')
  annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
  add <- rbind(add, annot)
  write.table(add, file = gct_file, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # phenotype file
  if(type == "between"){
    groups <- unique(meta$strain)
  } else {
    groups <- levels(factor(meta$label))
  }
  ngroups <- length(groups)
  ph <- matrix(nrow = 3, ncol = ncol(gct))
  # first row
  ph[1,1] <- ncol(gct)
  ph[1,2] <- ngroups
  ph[1,3] <- 1
  # second row
  ph[2,1:ngroups] <- groups
  ph[2,1] <- paste0('# ', ph[2,1])
  # third row
  if(type == "between"){
    ph[3,] <- meta$strain
  } else {
    ph[3,] <- meta$label
  }
  ph <- as.data.frame(ph)
  write.table(ph, file = cls_file, quote = F, sep = " ", na = "", col.names = F, row.names = F)
}

# across strains
gsea.input(counts_collapsed = expr.counts.mat,
           meta = meta_file,
           strains = unique(meta_file$strain),
           groups = c('exercised','non_exercised'),
           type = "across",
           gct_file = file.path(output_dir, paste0(prefix, '_gsea_all_strains.gct')),
           cls_file = file.path(output_dir, paste0(prefix, '_gsea_all_strains.cls')))

# within strain
strains <- unique(meta_file$strain)
for(i in 1:length(strains)){
  st <- strains[i]
  fname <- paste0('gsea_',st)
  gsea.input(counts_collapsed = expr.counts.mat,
             meta = meta_file,
             strains = st,
             groups = c('exercised','non_exercised'),
             type = "within",
             gct_file = file.path(output_dir, paste0(prefix, "_", fname, '.gct')),
             cls_file = file.path(output_dir, paste0(prefix, "_", fname, '.cls')))
}

# between strain 
# exercised
comb_strains <- as.data.frame(t(combn(strains, m = 2)))
colnames(comb_strains) <- c("strain1", "strain2")

for(i in 1:nrow(comb_strains)){
  st <- c(as.character(comb_strains[i,])) 
  fname <- paste('gsea', st[1], 'vs', st[2], 'exercised', sep = "_")
  gsea.input(counts_collapsed = expr.counts.mat,
             meta = meta_file,
             strains = st,
             groups = 'exercised',
             type = 'between',
             gct_file = file.path(output_dir, paste0(prefix, "_", fname, '.gct')),
             cls_file = file.path(output_dir, paste0(prefix, "_", fname, '.cls')))
}

# non_exercised
for(i in 1:nrow(comb_strains)){
  st <- c(as.character(comb_strains[i,])) 
  fname <- paste('gsea', st[1], 'vs', st[2], 'non_exercised', sep = "_")
  gsea.input(counts_collapsed = expr.counts.mat,
             meta = meta_file,
             strains = st,
             groups = 'non_exercised',
             type = 'between',
             gct_file = file.path(output_dir, paste0(prefix, "_", fname, '.gct')),
             cls_file = file.path(output_dir, paste0(prefix, "_", fname, '.cls')))
}
