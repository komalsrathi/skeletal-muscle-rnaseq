# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/mouse_skm/skm_mm10.RData \
--annot data/gencode.vM17.annotation.txt \
--prefix skm \
--outdir data/mouse_skm

# Batch correct for two batches (count matrix)
Rscript R/batch_correct.R \
--combined_mat data/mouse_skm/skm_collapsed_counts_matrix.RData \
--combined_clin data/mouse_skm/skm-meta-data.txt \
--corrected_outfile skm_collapsed_counts_matrix_batchcorrected.RData \
--type counts \
--outdir data/mouse_skm

# Batch correct for two batches (fpkm matrix)
Rscript R/batch_correct.R \
--combined_mat data/mouse_skm/skm_collapsed_fpkm_matrix.RData \
--combined_clin data/mouse_skm/skm-meta-data.txt \
--corrected_outfile skm_collapsed_fpkm_matrix_batchcorrected.RData \
--type fpkm \
--outdir data/mouse_skm

# Within strain exercised vs nonexercised
Rscript R/diffexpr_within_strain.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--type 'non_exercised' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix within-strain-ex-vs-nonex \
--excel TRUE \
--text TRUE \
--outdir results/skm/within-strain-ex-vs-nonex

# Between strain comparison of non-exercised mice
Rscript R/diffexpr_between_strain.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--type 'non_exercised' \
--col 'strain' \
--fc 0 \
--plx 7 \
--prefix between-strain-nonex \
--excel TRUE \
--text TRUE \
--outdir results/skm/between-strain-nonex

# Between strain comparison of exercised mice (responders + non-responders)
Rscript R/diffexpr_between_strain.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--type 'exercised' \
--col 'strain' \
--fc 0 \
--plx 7 \
--prefix between-strain-ex \
--excel TRUE \
--text TRUE \
--outdir results/skm/between-strain-ex

# All strain comparison of exercised mice (responders vs non-responders)
Rscript R/diffexpr_all_strains.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--type 'exercised_responders, exercised_non_responders' \
--col 'label' \
--fc 0 \
--plx 7 \
--prefix all-strains-exres-vs-exnonres \
--excel TRUE \
--text TRUE \
--outdir results/skm/all-strains-exres-vs-exnonres

# t-SNE sample clustering
Rscript R/clustering.R \
--fpkm_matrix data/mouse_skm/skm_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--outdir results/skm

# sample correlation
Rscript R/sample_correlation.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--outdir results/skm

# expression per sample (skm specific genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_skm/skm_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--gene_list data/gene-lists/skm-specific-genes.txt \
--title "Skeletal Muscle Specific Genes" \
--prefix skm_specific_genes_expr \
--group FALSE \
--view sample \
--width 15 \
--height 6 \
--outdir results/skm/expression-plots

# expression per strain (skm specific genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_skm/skm_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--gene_list data/gene-lists/skm-specific-genes.txt \
--title "Skeletal Muscle Specific Genes" \
--prefix skm_specific_genes_expr_byStrain \
--group FALSE \
--view strain \
--width 15 \
--height 6 \
--outdir results/skm/expression-plots

# expression per sample (housekeeping genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_skm/skm_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--gene_list data/gene-lists/house-keeping-genes.txt \
--title "House Keeping Genes" \
--prefix house_keeping_genes_expr \
--group FALSE \
--view sample \
--width 15 \
--height 6 \
--outdir results/skm/expression-plots

# expression per strain (housekeeping genes)
Rscript R/expression_boxplots.R \
--fpkm_matrix data/mouse_skm/skm_collapsed_fpkm_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--gene_list data/gene-lists/house-keeping-genes.txt \
--title "House Keeping Genes" \
--prefix house_keeping_genes_expr_byStrain \
--group FALSE \
--view strain \
--width 15 \
--height 6 \
--outdir results/skm/expression-plots

