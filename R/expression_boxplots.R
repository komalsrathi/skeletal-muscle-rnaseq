library(optparse)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(plyr)

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
  make_option(c("--group"), type = "character",
              help = "Use grouping var as facet? TRUE or FALSE"),
  make_option(c("--view"), type = "character",
              help = "x axis as sample or strain?"),
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
title <- opt$title
group <- opt$group
width <- as.numeric(opt$width)
height <- as.numeric(opt$height)
outdir <- opt$outdir
view <- opt$view

# create output directory
dir.create(outdir, recursive = TRUE, showWarnings = F)

# Expression plot of heart-specific genes
genes <- read.delim(gene_list, stringsAsFactors = F)
colnames(genes)[1] <- 'gene_symbol'
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

# add view
expr.fpkm <- expr.fpkm %>%
  dplyr::mutate('view' = !!as.name(view))
  
# generalized function
create.boxplot <- function(fpkm_mat, gene_list, title){
  
  print(unique(gene_list$type))
  
  if(length(unique(gene_list$type)) == 1 && unique(gene_list$type) == ""){
    print("breaking out of function...")
    return(NULL)
  }
  
  # subset to genes of interest
  fpkm_mat <- fpkm_mat %>%
    inner_join(gene_list, by = c('gene_symbol'))
  
  # set title
  if(is.null(title)){
    title <- unique(fpkm_mat$type)
  }
  
  # refactor
  syms <- fpkm_mat %>%
    group_by(view) %>%
    mutate(median = median(fpkm)) %>%
    select(view, median) %>%
    unique() %>%
    arrange(median) %>%
    .$view
  fpkm_mat$view <- factor(fpkm_mat$view, levels = syms)
  
  # convert fpkm to log2 fpkm
  fpkm_mat <- fpkm_mat %>%
    mutate(fpkm = log2(fpkm + 1))
  
  # boxplot
  p <- ggplot(fpkm_mat, aes(view, fpkm, color = strain)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
    ylab("log2 FPKM") + xlab(view) +
    theme_Publication2(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle(title) +
    stat_compare_means(method = "anova", label.y = max(fpkm_mat$fpkm) + 0.5,
                       label.x.npc = "center", label.y.npc = "top") +  
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")
  return(p)
}

if(group == TRUE){
  res <- dlply(.data = genes, 
               .variables = 'type', 
               .fun = function(x) create.boxplot(fpkm_mat = expr.fpkm, gene_list = x, title = NULL))
} else {
  res <- create.boxplot(fpkm_mat = expr.fpkm, gene_list = genes, title = title)
  res <- list(res)
}
res <- res[!sapply(res,is.null)]

# save output
fname <- file.path(outdir, paste0(prefix, '.pdf'))
glist <- lapply(res, ggplotGrob)
ggsave(marrangeGrob(grobs = glist, nrow=1, ncol=1), filename = fname, width = width, height = height)

