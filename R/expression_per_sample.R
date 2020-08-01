library(optparse)
library(ggplot2)
library(reshape2)
library(tidyverse)
source('~/Projects/Utils/pubTheme.R')

option_list <- list(
  make_option(c("--fpkm_matrix"), type = "character",
              help = "RData object of FPKM"),
  make_option(c("--meta_file"), type = "character",
              help = "Meta data"),
  make_option(c("--gene_list"), type = "character",
              help = "Gene list"),
  make_option(c("--title"), type = "character",
              help = "Title for output file"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output file"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
fpkm_matrix <- opt$fpkm_matrix
gene_list <- opt$gene_list
meta_file <- opt$meta_file
prefix <- opt$prefix
title <- opt$title
outdir <- opt$outdir

# Expression plot of heart-specific genes
genes <- read.delim(gene_list, header = F, stringsAsFactors = F)
genes <- unique(genes$V1)

# FPKM expression
load(fpkm_matrix)
expr.fpkm <- expr.fpkm[[1]]
expr.fpkm <- expr.fpkm[rownames(expr.fpkm) %in% genes,]
expr.fpkm <- melt(as.matrix(expr.fpkm))

# meta data
meta <- read.delim(meta_file, stringsAsFactors = F)
expr.fpkm <- expr.fpkm %>%
  inner_join(meta, by = c('Var2' = 'sample'))

# refactor
syms <- expr.fpkm %>%
  group_by(Var2) %>%
  mutate(median = median(value)) %>%
  select(Var2, median) %>%
  unique() %>%
  arrange(median) %>%
  .$Var2
expr.fpkm$Var2 <- factor(expr.fpkm$Var2, levels = syms)

# boxplot per sample
p <- ggplot(expr.fpkm, aes(Var2, log2(value + 1), color = strain)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  ylab("log2 FPKM") + xlab("Sample") +
  theme_Publication2(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle(title)
fname <- file.path(outdir, paste0(prefix,'.pdf'))
ggsave(p, filename = fname, width = 14, height = 6)

# boxplot per sample
# ggplot(expr.fpkm, aes(Var2, log2(value + 1))) + 
#   geom_boxplot()
