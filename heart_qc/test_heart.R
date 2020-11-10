setwd('~/Projects/Skeletal_Muscle_Patrick/')
library(tidyverse)
library(Rtsne)
source('R/filterExpr.R')

# metadata
dat <- read.delim('data/mouse_heart/heart_covariable_list.txt')
meta <- read.delim('data/mouse_heart/heart-meta-data.txt')
dat <- dat %>%
  inner_join(meta, by = c("mice" = "sample")) %>%
  arrange(strain) %>%
  mutate(tmp = mice) %>%
  column_to_rownames('tmp')

# expression
load('data/mouse_heart/heart_collapsed_counts_matrix.RData')
expr.counts = expr.counts[[1]]
expr.counts = expr.counts[,rownames(dat)]

# design
var <- factor(dat[,'strain'])
design <- model.matrix(~0+var)
colnames(design) <- levels(var)
rownames(design) <- dat$mice
print(dim(design))

# filter
identical(rownames(dat), colnames(expr.counts))
expr.counts = filterExpr(expr.counts.mat = expr.counts, group = dat[,'strain'], var.filter = F)

# normalize
y <- DGEList(counts = as.matrix(expr.counts))
y <- calcNormFactors(y)
v <- voom(counts = y, design = design, plot = FALSE)
voomData <- v$E

# clustering
set.seed(42)
tsne.out <- Rtsne(X = t(voomData), initial_dims = 2, perplexity = 3, max_iter = 1000)
tsne.out <- data.frame(tsne.out$Y, dat)
ggplot(tsne.out, aes(X1, X2)) +
  geom_point(size = 10, alpha = 0.5, aes(color = strain, shape = label)) +
  geom_text(aes(label = mice), size = 4) + 
  theme_bw() +
  ggtitle("Expressed genes: 15138 genes") + xlab("PC1") + ylab("PC2")

#  pca
prData <- prcomp(voomData)
pca.data <- prData$rotation
pca.data <- data.frame(pca.data)[1:4]
pca.data <- data.frame(pca.data, dat)
p = ggplot(pca.data, aes(PC1, PC2)) +
  geom_point(size = 5, alpha = 0.5, aes(color = strain, shape = label)) +
  geom_text(aes(label = mice), size = 2) + 
  theme_bw() +
  ggtitle("Expressed genes: 15138 genes") +
  geom_hline(yintercept = 0, lty = 2, col="grey50") +
  geom_vline(xintercept = 0, lty = 2, col="grey50") + 
  xlab(paste0("PC1: ",summary(prData)$importance[2,1]*100, "%")) +
  ylab(paste0("PC2: ",summary(prData)$importance[2,2]*100, "%"))

q = ggplot(pca.data, aes(PC2, PC3)) +
  geom_point(size = 5, alpha = 0.5, aes(color = strain, shape = label)) +
  geom_text(aes(label = mice), size = 2) + 
  theme_bw() +
  ggtitle("Expressed genes: 15138 genes") +
  geom_hline(yintercept = 0, lty = 2, col="grey50") +
  geom_vline(xintercept = 0, lty = 2, col="grey50") +
  xlab(paste0("PC2: ",summary(prData)$importance[2,2]*100, "%")) +
  ylab(paste0("PC3: ",summary(prData)$importance[2,3]*100, "%"))

library(ggpubr)
ggsave(filename = "heart_qc/pca_plot_all_strains.png", 
       plot = ggarrange(p, q, common.legend = T), device = "png", width = 15, height = 6)

# split into two groups
tmp = pca.data[order(pca.data$PC1, pca.data$PC2),]
tmp[1:35,'batch'] <- 'group1'
tmp$batch[tmp$batch != "group1"] = "group2"

# ggbiplot
# library(ggbiplot)
# prData <- prcomp(t(voomData))
# ggbiplot(prData, obs.scale = 1, var.scale = 1, 
#          groups = factor(tmp$batch),
#          ellipse = TRUE, circle = TRUE) +
#   scale_color_discrete(name = '') +
#   theme(legend.direction = 'horizontal', legend.position = 'top')
# 
# summary(prData)

# new design
var <- factor(tmp[,'batch'])
design <- model.matrix(~0+var)
colnames(design) <- levels(var)
rownames(design) <- tmp$mice
print(dim(design))
contrast.matrix <- makeContrasts(group1-group2, levels=design) 

# arrange expr
voomData <- voomData[,rownames(design)]

# now do differential expression
fit <- lmFit(voomData, design)
fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
fit2 <- eBayes(fit2)
outputLimma <- topTable(fit2, coef = 1, number = Inf)
outputLimma <- outputLimma %>%
  rownames_to_column('gene')
outputLimma$DEG  = ifelse(outputLimma$adj.P.Val < 0.05 & abs(outputLimma$logFC) > 1, "Yes", "No")
write.table(outputLimma, file = 'heart_qc/diffexpr_genes.txt', quote = F, sep = "\t", row.names = F)
outputLimma <- outputLimma %>%
  filter(DEG == "Yes") %>%
  dplyr::select(gene, logFC) %>%
  arrange(desc(logFC))

# use countdata and convert to human genes
expr.counts = expr.counts[,rownames(design)]

library(biomaRt)
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
human.genes <- convertMouseGeneList(rownames(expr.counts))
human.genes <- human.genes %>%
  filter(HGNC.symbol != "")
expr.counts.gsea <- expr.counts %>%
  rownames_to_column("gene_symbol") %>%
  inner_join(human.genes, by = c("gene_symbol" = "MGI.symbol")) %>%
  dplyr::select(-c(gene_symbol, Gene.type)) %>%
  dplyr::mutate(means = rowMeans(.[1:60])) %>%
  arrange(desc(means)) %>%
  distinct(HGNC.symbol, .keep_all = TRUE) %>% 
  dplyr::select(-c(means)) %>%
  column_to_rownames(var = "HGNC.symbol")
expr.counts.gsea = expr.counts.gsea[,rownames(design)]

# gsea with c8 data
library(GSEABase)
gene_set <- getGmt(file.path('heart_qc/c8.all.v7.2.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)

# try gage witih full expression
gene_set_out_group1 = gage(exprs = log2(expr.counts.gsea + 1), gsets = gene_set, ref = c(1:35), samp = c(36:60), compare = "unpaired")
gene_set_out_group1.up = gene_set_out_group1$greater %>% as.data.frame %>% filter(q.val < 0.05)
gene_set_out_group1.down = gene_set_out_group1$less %>% as.data.frame %>% filter(q.val < 0.05)
# gene_set_out_group2 = gage(exprs = log2(expr.counts.gsea + 1), gsets = gene_set, samp = c(1:35), ref = c(36:60), compare = "unpaired")
# gene_set_out_group2.up = gene_set_out_group2$greater %>% as.data.frame %>% filter(q.val < 0.05)
# gene_set_out_group2.down = gene_set_out_group2$less %>% as.data.frame %>% filter(q.val < 0.05)
gene_set_out_group1.down = gene_set_out_group1.down  %>% rownames_to_column('pathway') %>% dplyr::select(pathway, p.val, set.size)
gene_set_out_group1.up = gene_set_out_group1.up  %>% rownames_to_column('pathway') %>% dplyr::select(pathway, p.val, set.size)
write.table(gene_set_out_group1.down, file = 'heart_qc/gage_gene_set_out_group1_down.txt', quote = F, sep = "\t", row.names = F)
write.table(gene_set_out_group1.up, file = 'heart_qc/gage_gene_set_out_group1_up.txt', quote = F, sep = "\t", row.names = F)

# try fgsea on preranked list of group a vs group b
library(fgsea)
outputLimma.human = human.genes %>%
  inner_join(outputLimma, by = c("MGI.symbol" = "gene"))
gene_list = outputLimma.human$logFC
names(gene_list) = outputLimma.human$HGNC.symbol
gene_list = sort(gene_list, decreasing = TRUE)
myGO = fgsea::gmtPathways('heart_qc/c8.all.v7.2.symbols.gmt')
fgRes <- fgsea(pathways = myGO, 
               stats = gene_list,
               minSize = 15,
               maxSize = 1500,
               nperm=10000) %>% 
  as.data.frame() %>%
  filter(pval < 0.05)
fgRes = fgRes %>%
  mutate(genes = sapply(leadingEdge, paste, collapse=",")) %>%
  dplyr::select(-c(leadingEdge)) %>%
  mutate(Direction = ifelse(NES < 0, "Down", "Up"))
write.table(fgRes, file = 'heart_qc/fgsea_ranked_pathway_enrichment.txt', quote = F, sep = "\t", row.names = F)

# try gsva with full expression
library(GSVA)
ssgsea_scores <- GSVA::gsva(as.matrix(log2(expr.counts.gsea + 1)),
                          gene_set,
                          method = "ssgsea",
                          min.sz=1, max.sz=1500, 
                          mx.diff = TRUE)

# plot these by combining groups
library(reshape2)
ssgsea_scores <- melt(ssgsea_scores, varnames = c("pathway","mice"), value.name = "ssgsea_score")
ssgsea_scores <- ssgsea_scores %>%
  inner_join(tmp, by = "mice")
ssgsea_scores %>% group_by(pathway) %>% mutate(test = t.test(ssgsea_score))

pathways = unique(ssgsea_scores$pathway)
vec = ''
for(i in 1:length(pathways)){
  print(i)
  ssgsea_scores_tmp = ssgsea_scores %>%
    filter(pathway == pathways[i])
  one = ssgsea_scores_tmp %>%
    filter(batch == "group1") %>%
    .$ssgsea_score
  two = ssgsea_scores_tmp %>%
    filter(batch == "group2") %>%
    .$ssgsea_score
  t_test_out = t.test(one, two, paired = FALSE)
  pval = t_test_out$p.value
  print(pathways[i])
  print(pval)
  if(pval < 0.01){
    vec = c(vec, pathways[i])
  }
}

ggplot(ssgsea_scores, aes(x = pathway, y = ssgsea_score, color = batch)) +
  facet_wrap(~pathway) +
  geom_boxplot() +
  stat_compare_means(method = "t.test")
