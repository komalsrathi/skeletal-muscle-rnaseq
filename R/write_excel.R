# send gene expression on the individual mouse level
# mouse heart
load('data/mouse_heart/heart_collapsed_fpkm_matrix.RData')
expr.fpkm <- expr.fpkm[[1]]
load('data/mouse_heart/heart_collapsed_counts_matrix.RData')
expr.counts <- expr.counts[[1]]
writexl::write_xlsx(x = data.frame(gene = rownames(expr.counts), expr.counts, check.names = F), path = 'data/mouse_heart/heart_collapsed_counts_matrix.xlsx', col_names = TRUE, format_headers = TRUE)
writexl::write_xlsx(x = data.frame(gene = rownames(expr.fpkm), expr.fpkm, check.names = F), path = 'data/mouse_heart/heart_collapsed_fpkm_matrix.xlsx', col_names = TRUE, format_headers = TRUE)

# mouse skm
load('data/mouse_skm/skm_collapsed_fpkm_matrix.RData')
expr.fpkm <- expr.fpkm[[1]]
load('data/mouse_skm/skm_collapsed_counts_matrix.RData')
expr.counts <- expr.counts[[1]]
writexl::write_xlsx(x = data.frame(gene = rownames(expr.counts), expr.counts, check.names = F), path = 'data/mouse_skm/skm_collapsed_counts_matrix.xlsx', col_names = TRUE, format_headers = TRUE)
writexl::write_xlsx(x = data.frame(gene = rownames(expr.fpkm), expr.fpkm, check.names = F), path = 'data/mouse_skm/skm_collapsed_fpkm_matrix.xlsx', col_names = TRUE, format_headers = TRUE)
