# Author: Komal S. Rathi
# Function: PCA plot of voom normalized data
# Date: 04/03/2020

# Generalized function
source('../Utils/pubTheme.R')

library(Rtsne)
library(ggpubr)

tsne.plot <- function(voomData, meta, fname, plx){
  
  # subset voom normalized data using meta file
  rownames(meta) <- meta$sample
  voomData <- voomData[,rownames(meta)]
  
  # check if meta and expression are compatible
  if(identical(rownames(meta), colnames(voomData))) {
    print("Proceed")
  } else {
    break
  }
  
  # t-SNE before combat adjustment
  set.seed(42)
  tsneOut <- Rtsne(t(voomData), initial_dims = 50, perplexity = plx, max_iter = 1000)
  tsneOut <- data.frame(tsneOut$Y, meta)
  p <- ggplot(tsneOut, aes(X1, X2)) +
    geom_point(size = 5, alpha = 0.5, aes(color = strain, shape = label)) +
    geom_text(aes(label = sample), size = 2) + 
    theme_bw() +
    ggtitle("T-SNE Clustering (Voom normalized data)") +
    theme_Publication2() + xlab("PC1") + ylab("PC2") 
  ggsave(filename = fname, plot = p, device = "pdf", width = 7, height = 5)
}

pca.plot <- function(voomData, meta, fname){
  
  # subset voom normalized data using meta file
  rownames(meta) <- meta$sample
  voomData <- voomData[,rownames(meta)]
  
  # check if meta and expression are compatible
  if(identical(rownames(meta), colnames(voomData))) {
    print("Proceed")
  } else {
    break
  }
  
  # t-SNE before combat adjustment
  prData <- prcomp(voomData)
  pca.data <- prData$rotation
  pca.data <- data.frame(pca.data)[1:4]
  pca.data <- data.frame(pca.data, meta)
  pdf(file = fname, width = 10, height = 5, onefile = FALSE)
  p <- ggplot(pca.data, aes(PC1, PC2)) +
    geom_point(size = 5, alpha = 0.5, aes(color = strain, shape = label)) +
    geom_text(aes(label = sample), size = 2) + 
    theme_bw() +
    ggtitle("PCA Clustering (Voom normalized data)") +
    theme_Publication2()  
  q <- ggplot(pca.data, aes(PC3, PC4)) +
    geom_point(size = 5, alpha = 0.5, aes(color = strain, shape = label)) +
    geom_text(aes(label = sample), size = 2) + 
    theme_bw() +
    ggtitle("PCA Clustering (Voom normalized data)") +
    theme_Publication2() 
  print(ggarrange(p, q, common.legend = T))
  dev.off()
}
