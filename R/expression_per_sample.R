library(optparse)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
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
  make_option(c("--sample"), type = "character",
              help = "x-axis by sample? TRUE or FALSE"),
  make_option(c("--strain"), type = "character",
              help = "x-axis by strain? TRUE or FALSE"),
  make_option(c("--group"), type = "character",
              help = "Use grouping var as facet? TRUE or FALSE"),
  make_option(c("--view"), type = "character",
              help = "x axis as sample or strain? TRUE or FALSE"),
  make_option(c("--width"), type = "character",
              help = "Width of plot"),
  make_option(c("--height"), type = "character",
              help = "Height of plot"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
fpkm_matrix <- opt$fpkm_matrix
gene_list <- opt$gene_list
meta_file <- opt$meta_file
prefix <- opt$prefix
sample <- opt$sample
strain <- opt$strain
title <- opt$title
group <- opt$group
w <- as.numeric(opt$width)
h <- as.numeric(opt$height)
outdir <- opt$outdir
view <- opt$view

# Expression plot of heart-specific genes
genes <- read.delim(gene_list, stringsAsFactors = F)
genes <- unique(genes)
gene.sym <- unique(genes$gene_symbol)

# FPKM expression
load(fpkm_matrix)
expr.fpkm <- expr.fpkm[[1]]
expr.fpkm <- expr.fpkm[rownames(expr.fpkm) %in% gene.sym,]
expr.fpkm <- melt(as.matrix(expr.fpkm), varnames = c("gene_symbol", "sample"), value.name = 'fpkm')

# meta data
meta <- read.delim(meta_file, stringsAsFactors = F)
expr.fpkm <- expr.fpkm %>%
  inner_join(meta, by = c('sample'))

# if group is TRUE, plot cell type groups
if(group == TRUE){
  
  # merge grouping of genes
  expr.fpkm <- expr.fpkm %>%
    inner_join(genes, by = c('gene_symbol')) %>%
    filter(cell_type != "")
  
  celltypes <- unique(expr.fpkm$cell_type)
  newList <- list()
  for(i in 1:length(celltypes)){
    dat <- expr.fpkm %>%
      filter(cell_type %in% celltypes[i]) %>%
      mutate(fpkm = log2(fpkm + 1))
    p <- ggplot(dat, aes_string(x = view, y = 'fpkm', color = 'strain')) + 
      stat_boxplot(geom ='errorbar', width = 0.2) +
      geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
      ylab("log2 FPKM") + xlab("Sample") +
      theme_Publication2(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      ggtitle(celltypes[i]) +
      stat_compare_means(method = "anova", label.y = max(dat$fpkm)+0.5,
                         label.x.npc = "center", label.y.npc = "top") +  
      stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")
    newList[[i]] <- p
  }
}

# boxplot per sample
if(sample == TRUE){
  
  # refactor
  syms <- expr.fpkm %>%
    group_by(sample) %>%
    mutate(median = median(fpkm)) %>%
    select(sample, median) %>%
    unique() %>%
    arrange(median) %>%
    .$sample
  expr.fpkm$sample <- factor(expr.fpkm$sample, levels = syms)
  newList <- list()
  p <- ggplot(expr.fpkm, aes(sample, log2(fpkm + 1), color = strain)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
    ylab("log2 FPKM") + xlab("Sample") +
    theme_Publication2(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title) +
    stat_compare_means(method = "anova", label.y = max(log2(expr.fpkm$fpkm + 1))+0.5,
                       label.x.npc = "center", label.y.npc = "top") +  
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")
  newList[[1]] <- p
} 

# boxplot per strain
if(strain == TRUE){
  
  # refactor
  syms <- expr.fpkm %>%
    group_by(strain) %>%
    mutate(median = median(fpkm)) %>%
    select(strain, median) %>%
    unique() %>%
    arrange(median) %>%
    .$strain
  expr.fpkm$strain <- factor(expr.fpkm$strain, levels = syms)
  newList <- list()
  p <- ggplot(expr.fpkm, aes(strain, log2(fpkm + 1), color = strain)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
    ylab("log2 FPKM") + xlab("Sample") +
    theme_Publication2(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title) +
    stat_compare_means(method = "anova", label.y = max(log2(expr.fpkm$fpkm + 1))+0.5,
                       label.x.npc = "center", label.y.npc = "top") +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")
  newList[[1]] <- p
}

fname <- file.path(outdir, paste0(prefix,'.pdf'))
glist <- lapply(newList, ggplotGrob)
ggsave(marrangeGrob(grobs = glist, nrow=1, ncol=1), filename = fname, width = w, height = h)

