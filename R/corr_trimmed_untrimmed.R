# Author: Komal S. Rathi
# Date: 03/09/2020
# Function: Correlation of FPKM expression in trimmed and untrimmed data

setwd('~/Projects/Skeletal_Muscle_Patrick/')

A13 <- data.table::fread('data/A13.genes.results')
A13 <- A13[,c('gene_id', 'FPKM')]
colnames(A13)[2] <- 'trimmed'
load('data/skeletal_muscle_mm10.RData')
expr.fpkm <- expr.fpkm[,c('gene_id','A13')]
colnames(expr.fpkm)[2] <- 'orig'
total <- merge(A13, expr.fpkm, by = 'gene_id')

load('data/mm10_annotation.RData')
total <- merge(total, annot, by = 'gene_id')
total <- total[grep('^Mt', total$gene_symbol),]

cor(total$trimmed, total$orig, method = "pearson") # 0.99
cor(total$trimmed, total$orig, method = "spearman") # 0.99
plot(total$trimmed, total$orig)
