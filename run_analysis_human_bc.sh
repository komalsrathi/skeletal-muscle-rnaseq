# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/human/human_hg38.RData \
--annot data/gencode.v23.annotation.txt \
--prefix human \
--outdir data/human

# Batch correct for two batches (count matrix)
Rscript R/batch_correct.R \
--combined_mat data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--combined_clin data/human/human-meta-data.txt \
--corrected_outfile human_collapsed_counts_matrix_batchcorrected.RData \
--type counts \
--outdir data/human

# Batch correct for two batches (fpkm matrix)
Rscript R/batch_correct.R \
--combined_mat data/human/human_collapsed_fpkm_matrix.RData \
--combined_clin data/human/human-meta-data.txt \
--corrected_outfile human_collapsed_fpkm_matrix_batchcorrected.RData \
--type fpkm \
--outdir data/human

## Variance filter
# 399_2 vs 370_1
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter TRUE \
--type '399_2, 370_1' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_1' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr-399_2-vs-370_1

# 399_2 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter TRUE \
--type '399_2, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr-399_2-vs-370_5

# 370_1 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter TRUE \
--type '370_1, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-370_1-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr-370_1-vs-370_5

## Remove filter by variance
# 399_2 vs 370-1
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter FALSE \
--type '399_2, 370_1' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_1' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr-399_2-vs-370_1_novarfilter

# 399_2 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter FALSE \
--type '399_2, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr-399_2-vs-370_5_novarfilter

# 370_1 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter FALSE \
--type '370_1, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-370_1-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr-370_1-vs-370_5_novarfilter

# remove variance filter and do one comparison tumor vs normals
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter FALSE \
--type 'Treat, Control' \
--col 'treat' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr_treat_vs_control' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr_treat_vs_control_novarfilter

# with variance filter
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix_batchcorrected.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/human/human-meta-data.txt \
--var_filter TRUE \
--type 'Treat, Control' \
--col 'treat' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr_treat_vs_control' \
--excel TRUE \
--text TRUE \
--outdir results/human_bc/diffexpr_treat_vs_control