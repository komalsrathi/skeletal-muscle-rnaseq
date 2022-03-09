library(tidyverse)

# skeletal muscle
# counts
load('data/mouse_skm/skm_collapsed_counts_matrix.RData')
expr_counts <- expr.counts[[1]]
colnames(expr_counts)
dat <- read.delim('data/mouse_skm/skm-meta-data.txt')
dat$sample <- as.character(dat$sample)
expr_counts <- expr_counts %>%
  dplyr::select(dat$sample)
expr_counts <- expr_counts %>%
  rownames_to_column("Gene_Symbol")
readr::write_tsv(expr_counts, file = 'data/mouse_skm/skeletal_muscle_expected_counts_matrix.tsv')

# fpkm
load('data/mouse_skm/skm_collapsed_fpkm_matrix.RData')
expr_fpkm <- expr.fpkm[[1]]
colnames(expr_fpkm)
dat <- read.delim('data/mouse_skm/skm-meta-data.txt')
dat$sample <- as.character(dat$sample)
expr_fpkm <- expr_fpkm %>%
  dplyr::select(dat$sample)
expr_fpkm <- expr_fpkm %>%
  rownames_to_column("Gene_Symbol")
readr::write_tsv(expr_fpkm, file = 'data/mouse_skm/skeletal_muscle_fpkm_matrix.tsv')

# Heart
# counts
load('data/mouse_heart/heart_collapsed_counts_matrix.RData')
expr_counts <- expr.counts[[1]]
colnames(expr_counts)
dat <- read.delim('data/mouse_heart/heart-meta-data.txt')
dat$sample <- as.character(dat$sample)
expr_counts <- expr_counts %>%
  dplyr::select(dat$sample)
expr_counts <- expr_counts %>%
  rownames_to_column("Gene_Symbol")
readr::write_tsv(expr_counts, file = 'data/mouse_heart/heart_expected_counts_matrix.tsv')

# fpkm
load('data/mouse_heart/heart_collapsed_fpkm_matrix.RData')
expr_fpkm <- expr.fpkm[[1]]
colnames(expr_fpkm)
dat <- read.delim('data/mouse_heart/heart-meta-data.txt')
dat$sample <- as.character(dat$sample)
expr_fpkm <- expr_fpkm %>%
  dplyr::select(dat$sample)
expr_fpkm <- expr_fpkm %>%
  rownames_to_column("Gene_Symbol")
readr::write_tsv(expr_fpkm, file = 'data/mouse_heart/heart_fpkm_matrix.tsv')
