# Author: Komal S. Rathi
# Date: 03/05/2020
# Function: % Unique Mapped Reads

library(ggplot2)
source('~/Projects/Utils/pubTheme.R')

perc.aligned <- read.delim('data/aligned_perc.txt', stringsAsFactors = F)

df <- read.delim('data/meta-data.txt', stringsAsFactors = F)
df$type <- gsub("[0-9].*", "", df$sample)
rownames(df) <- df$sample
df <- merge(perc.aligned, df, by = 'sample')

df$percent_aligned <- gsub("%", "", df$percent_aligned)
df$percent_aligned <- as.numeric(df$percent_aligned)
df <- df[order(df$percent_aligned),]
df$sample <- factor(df$sample, levels = df$sample)
p <- ggplot(df, aes(x = type, y = percent_aligned, fill = type)) +
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
  theme_Publication2(base_size = 12) +
  xlab("")  + ylab("% Unique Mapped Reads")
p
ggsave(p, filename = 'results/percent-aligned.pdf', width = 8, height = 6)

  
