# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/mouse_heart/heart_mm10.RData \
--annot data/gencode.vM17.annotation.txt \
--prefix heart \
--outdir data/mouse_heart

# Batch correct for two batches (count matrix)
Rscript R/batch_correct.R \
--combined_mat data/mouse_heart/heart_collapsed_counts_matrix.RData \
--combined_clin data/mouse_heart/heart-meta-data.txt \
--corrected_outfile heart_collapsed_counts_matrix_batchcorrected.RData \
--type counts \
--outdir data/mouse_heart

# Batch correct for two batches (fpkm matrix)
Rscript R/batch_correct.R \
--combined_mat data/mouse_heart/heart_collapsed_fpkm_matrix.RData \
--combined_clin data/mouse_heart/heart-meta-data.txt \
--corrected_outfile heart_collapsed_fpkm_matrix_batchcorrected.RData \
--type fpkm \
--outdir data/mouse_heart

# Within strain exercised vs nonexercised
Rscript R/diffexpr_within_strain.R \
--counts_matrix data/mouse_heart/heart_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--type 'non_exercised' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix within-strain-ex-vs-nonex \
--excel TRUE \
--text TRUE \
--outdir results/heart/within-strain-ex-vs-nonex

# Between strain comparison of non-exercised mice
Rscript R/diffexpr_between_strain.R \
--counts_matrix data/mouse_heart/heart_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--type 'non_exercised' \
--col 'strain' \
--fc 0 \
--plx 7 \
--prefix between-strain-nonex \
--excel TRUE \
--text TRUE \
--outdir results/heart/between-strain-nonex

# Between strain comparison of exercised mice (responders + non-responders)
Rscript R/diffexpr_between_strain.R \
--counts_matrix data/mouse_heart/heart_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--type 'exercised' \
--col 'strain' \
--fc 0 \
--plx 7 \
--prefix between-strain-ex \
--excel TRUE \
--text TRUE \
--outdir results/heart/between-strain-ex

# t-SNE sample clustering
Rscript R/clustering.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--outdir results/heart

# sample correlation
Rscript R/sample_correlation.R \
--counts_matrix data/mouse_heart/heart_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--outdir results/heart

# expression per sample (heart specific genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--gene_list data/gene-lists/heart-specific-genes.txt \
--title "Heart Specific Genes" \
--prefix heart_specific_genes_expr \
--group FALSE \
--view sample \
--width 15 \
--height 6 \
--outdir results/heart/expression-plots

# expression per strain (heart specific genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--gene_list data/gene-lists/heart-specific-genes.txt \
--title "Heart Specific Genes" \
--prefix heart_specific_genes_expr_byStrain \
--group FALSE \
--view strain \
--width 15 \
--height 6 \
--outdir results/heart/expression-plots

# # expression per sample (housekeeping genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--gene_list data/gene-lists/house-keeping-genes.txt \
--title "House Keeping Genes" \
--prefix house_keeping_genes_expr \
--group FALSE \
--view sample \
--width 15 \
--height 6 \
--outdir results/heart/expression-plots

# expression per strain (housekeeping genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--gene_list data/gene-lists/house-keeping-genes.txt \
--title "House Keeping Genes" \
--prefix house_keeping_genes_expr_byStrain \
--group FALSE \
--view strain \
--width 15 \
--height 6 \
--outdir results/heart/expression-plots

# expression per sample grouped by cell type per strain
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--gene_list data/gene-lists/heart-specific-genes.txt \
--title NULL \
--prefix heart_specific_genes_expr_byCelltypes_byStrain \
--group TRUE \
--view strain \
--width 15 \
--height 6 \
--outdir results/heart/expression-plots

# expression per sample grouped by cell type per sample
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_heart/heart_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--gene_list data/gene-lists/heart-specific-genes.txt \
--title NULL \
--prefix heart_specific_genes_expr_byCelltypes \
--group TRUE \
--view sample \
--width 15 \
--height 6 \
--outdir results/heart/expression-plots

# # Within strain ex vs nonex
# # using batch corrected data obtained with atrial cardiomyocyte expression levels
# Rscript R/diffexpr_within_strain.R \
# --counts_matrix data/mouse_heart/heart_collapsed_counts_matrix_batchcorrected.RData \
# --meta_file data/mouse_heart/heart-meta-data.txt \
# --type non_exercised \
# --prefix within-strain-ex-vs-nonex-batchcorrected \
# --outdir results/heart
