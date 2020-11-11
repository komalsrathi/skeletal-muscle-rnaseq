# Author: Komal S. Rathi
# Function: Script to generate input lists for GSEA preranked
setwd('~/Projects/Skeletal_Muscle_Patrick/')

# using GSEA preranked: 
# https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results

library(xlsx)
library(tidyverse)
library(biomaRt)

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

lf <- list.files('results/skeletal_muscle/covars/', pattern = ".xlsx", full.names = T)
lf <- grep('ANT1|IAI|EC77|B6_ME|weight|activity.CLAMS|VO2.end.max|VO2.avg.CLAMS|covar_association', lf, invert = T, value = T)

for(i in 1:length(lf)){
  dat <- read.xlsx(file = lf[i], sheetIndex = 1)
  dat <- dat %>%
    dplyr::select(gene, estimate, model.rSquared, p.value)
  
  # convert to human genes
  human.genes <- convertMouseGeneList(dat$gene)
  human.genes <- human.genes %>%
    filter(HGNC.symbol != "")
  res <- dat %>%
    inner_join(human.genes, by = c("gene" = "MGI.symbol")) %>%
    distinct(HGNC.symbol, .keep_all = TRUE) %>%
    dplyr::select(HGNC.symbol, p.value, estimate, model.rSquared)
  
  # version 1
  res.v1 <- res %>%
    arrange(p.value, desc(model.rSquared)) %>%
    mutate(rank = seq(from = 1, to = nrow(res), by = 1))
  
  # version 2
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
  res.v2 <- rbind(pos, neg)
  res.v2 <- res.v2 %>%
    arrange(rank)
  fname <- gsub('.*/|_model_output.xlsx', '', lf[i])
  fname <- file.path('results/gsea_preranked_input', fname)
  write.table(res.v1[,c('HGNC.symbol','rank')], file = paste0(fname,'_rankedlist_v1.txt'), quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(res.v2[,c('HGNC.symbol','rank')], file = paste0(fname,'_rankedlist_v2.txt'), quote = F, sep = "\t", row.names = F, col.names = F)
}

