# Author: Komal S. Rathi
# Date: 03/05/2020
# Function: Expression analysis of Mitochondrial Genes

library(optparse)
library(tidyverse)
library(reshape2)
source('~/Projects/Utils/pubTheme.R')

option_list <- list(
  make_option(c("--fpkm_matrix"), type = "character",
              help = "RData object of FPKM"),
  make_option(c("--gene_list"), type = "character",
              help = "Gene list"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
fpkm_matrix <- opt$fpkm_matrix
outdir <- opt$outdir

# load data
load(fpkm_matrix)
expr.fpkm <- expr.fpkm[[1]]

# house keeping gense
house.keeping.genes <- data.frame(gene_symbol = c("RPS18", "GAPDH", "PGK1", "PPIA", "RPL13A",
                                                  "RPLP0", "B2M", "YWHAZ", "SDHA", "TFRC",
                                                  "HPRT1"), type = "House-keeping")
house.keeping.genes$gene_symbol <- stringi::stri_trans_totitle(house.keeping.genes$gene_symbol)

# mito genes
mito.genes <- expr.fpkm[grep('^mt-', rownames(expr.fpkm)),]
mito.genes <- mito.genes %>%
  rownames_to_column('gene_symbol') %>%
  mutate(type = 'Mitochondria') %>%
  dplyr::select(gene_symbol, type)

# skeletal muscle genes
sm.genes <- data.frame(gene_symbol = c("Tnnt3", "Pvalb", "Mylpf", "Tnnc2", "Tnni2", "Myl1", "Tpm1", "Atp2a1", "Myh4", "Acta1"),
                       type = "Skeletal Muscle Specific")

# total genes
genes <- rbind(house.keeping.genes, mito.genes, sm.genes)
genes.expr <- expr.fpkm %>%
  rownames_to_column('gene_symbol') %>%
  inner_join(genes, by = c("gene_symbol"))
genes.expr <- melt(genes.expr)

# 15/37 expressed
genes.expr %>%
  group_by(gene_symbol) %>%
  mutate(median = median(value)) %>%
  select(gene_symbol, median, type) %>%
  filter(type == "Mitochondria") %>%
  filter(median > 0) %>%
  unique()

# refactor
syms <- genes.expr %>%
  group_by(gene_symbol) %>%
  mutate(median = median(value)) %>%
  select(gene_symbol, median) %>%
  unique() %>%
  arrange(median) %>%
  .$gene_symbol
genes.expr$gene_symbol <- factor(genes.expr$gene_symbol, levels = syms)
p <- ggplot(genes.expr, aes(gene_symbol, log2(value + 1), color = type)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  ylab("log2 FPKM") + xlab("Genes") +
  theme_Publication2(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Gene Expression")
p
fname <- file.path(outdir, 'FPKM_expression.pdf')
ggsave(p, filename = fname, width = 14, height = 6)

# barplot of CO1 genes in all E samples
co1 <- expr.fpkm %>%
  rownames_to_column('gene_symbol') %>%
  filter(gene_symbol == "mt-Co1") %>%
  melt()
co1$type <- gsub('[0-9]*', '',co1$variable)
p <- ggplot(co1, aes(variable, value, fill = type)) + 
  geom_bar(stat = "identity") + 
  ylab("FPKM") + xlab("Samples") +
  theme_Publication2(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("mt-Co1 Gene Expression")
q  <- ggplot(co1, aes(type, value, fill = type)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  ylab("FPKM") + xlab("Type") +
  theme_Publication2(base_size = 14) +
  ggtitle("mt-Co1 Gene Expression")
fname <- file.path(outdir, 'mtCo1_expression.pdf')
ggsave(gridExtra::grid.arrange(p, q, nrow = 2), filename = fname, width = 14, height = 14)

# intra group variability in CO1, CO2 and CO3 genes vs Skeletal muscle genes
# variability in mito genes vs house keeping and skeletal muscle
myCV <- function(x){ 
  sd(x)/mean(x)
}
allCVs <- expr.fpkm[apply(expr.fpkm, 1, function(x) !all(x==0)),]
allCVs <- apply(allCVs, FUN=myCV, MARGIN=1)
allCVs <- data.frame("CV" = allCVs)
allCVs <- allCVs %>%
  rownames_to_column("gene_symbol") %>%
  inner_join(genes.expr, by = c("gene_symbol"))
p <- ggplot(allCVs, aes(x = CV, fill = type)) + geom_density(alpha = 0.5) +
  theme_Publication(base_size = 12) +
  xlab("Coefficient of Variation") + ylab("Frequency")
fname <- file.path(outdir, 'CV_per_type.pdf')
ggsave(plot = p, filename = fname, width = 10, height = 6)

