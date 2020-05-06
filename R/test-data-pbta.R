# Author: Komal S. Rathi
# Date: 03/10/2020
# Function: Compare this data to PBTA mito expression
setwd('~/Projects/Skeletal_Muscle_Patrick/')

library(reshape2)
library(ggplot2)
source('~/Projects/Utils/pubTheme.R')

polya <- readRDS('~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm.polya.rds')
polya <- polya[grep('_MT-', polya$gene_id),]
polya$gene_id <- gsub('.*_','',polya$gene_id)
polya <- melt(polya)
polya$lab <- 'polyA'

stranded <- readRDS('~/Projects/OpenPBTA-analysis/data/pbta-gene-expression-rsem-fpkm.stranded.rds')
stranded <- stranded[grep('_MT-', stranded$gene_id),]
stranded$gene_id <- gsub('.*_','',stranded$gene_id)
stranded <- melt(stranded)
stranded$lab <- 'stranded'

total <- rbind(polya, stranded)
p <- ggplot(total, aes(x = gene_id, y = log2(value + 1))) +
  geom_boxplot() +
  facet_wrap(~lab, nrow = 2) +
  ggtitle("Test Dataset Comparison (Human)") +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust = 1))
p
ggsave(plot = p, filename = '~/Projects/Skeletal_Muscle_Patrick/results/PBTA_test_mitogenes.pdf', width = 10, height = 10)
