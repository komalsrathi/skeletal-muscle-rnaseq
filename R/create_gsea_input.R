# Author: Komal S. Rathi
# Function: Create GSEA input files
# Date: 05/26/2020

# libraries
library(tidyverse)
library(reshape2)
library(biomaRt)

# load count matrix
load('data/collapsed_counts_matrix.RData')
expr.counts.mat <- expr.counts[[1]]
expr.counts.mat <- expr.counts.mat[apply(expr.counts.mat!=0, 1, all),]

#expr.counts.annot <-  expr.counts[[2]]
meta <- read.delim('data/meta-data-withlitter.txt', stringsAsFactors = F)
meta$label <- gsub('-','_',meta$label)
meta$litter[meta$litter == ""] <- "U"
rownames(meta) <- meta$sample
meta$label2 <- ifelse(meta$label == "non_exercised", meta$label, "exercised")

# function to convert mouse to human
# function to convert data to human genes
# also output gene type (to do)
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
  dplyr::mutate(means = rowMeans(.[1:48])) %>%
  arrange(desc(means)) %>%
  distinct(HGNC.symbol, .keep_all = TRUE) %>% 
  dplyr::select(-c(means)) %>%
  column_to_rownames(var = "HGNC.symbol")

# function
# 1. Overall Ex vs Nonex: So only 2 groups
# 2. By strain Ex vs Nonex: So 8 groups
system('mkdir -p results/gsea/norm')

gsea.input <- function(counts.collapsed, meta, groupby, gct.file, cls.file) {
  # make file GSEA web version
  # expression file  -  counts file 
  # df <- data.frame(snames = colnames(counts.collapsed), colsplit(colnames(counts.collapsed), pattern = "_", names = c("sample","type")))
  if(groupby == "strain"){
    meta$strain <- gsub(' ','_',meta$strain)
    meta$label2 <- paste0(meta$strain, '_',meta$label2)
  }
  
  # order
  meta <- meta[order(meta$label2),]
  counts.collapsed <- counts.collapsed[,rownames(meta)]
  
  gct <- counts.collapsed[,rownames(meta)]
  add <- data.frame(NAME = c("#1.2",nrow(gct),"NAME"), Description = c('', ncol(gct), "Description"))
  total.cols <- ncol(gct) + 2
  add[,3:total.cols] <- ''
  colnames(add)[3:total.cols] <- colnames(gct)
  add[3,3:total.cols] <- colnames(gct)
  annot <- data.frame(NAME = rownames(gct), Description = 'na')
  annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
  add <- rbind(add, annot)
  write.table(add, file = gct.file, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # phenotype file
  groups <- levels(factor(meta$label2))
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
  ph[3,] <- meta$label2
  ph <- as.data.frame(ph)
  write.table(ph, file = cls.file, quote = F, sep = " ", na = "", col.names = F, row.names = F)
  
  # chp file
  # chp <- data.frame('Probe Set ID' = rownames(gct), 'Gene Symbol' = rownames(gct), 'Gene Title' = 'na', check.names = F)
  # write.table(chp, file = 'data/gsea_input/mapping.chip', quote = F, sep = "\t", row.names = F)
}

# by group
gsea.input(counts.collapsed = expr.counts.mat,
           meta = meta,
           groupby = 'label2',
           gct.file = 'results/gsea/matrix_byGroup.gct',
           cls.file = 'results/gsea/phenotype/phenotype_byGroup.cls')

# by strain
strains <- unique(meta$strain)
for(i in 1:length(strains)){
  st <- strains[i]
  strain.meta <- meta %>%
    rownames_to_column("id") %>%
    filter(strain == st) %>%
    column_to_rownames("id")
  strain.counts <- expr.counts.mat[,rownames(strain.meta)]
  gct.file <- paste0('results/gsea/matrix_byStrain_', st,'.gct')
  cls.file <- paste0('results/gsea/phenotype/phenotype_byStrain_', st,'.cls')
  gsea.input(counts.collapsed = strain.counts,
             meta = strain.meta,
             groupby = 'label2',
             gct.file = gct.file,
             cls.file = cls.file)
}
# gsea.input(counts.collapsed = expr.counts.mat,
#            meta = meta,
#            groupby = 'strain',
#            gct.file = 'results/gsea/matrix_byStrain.gct',
#            cls.file = 'results/gsea/phenotype_byStrain.cls')
