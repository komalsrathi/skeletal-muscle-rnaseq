# Merge RSEM
Rscript R/merge_rsem.R \
--sourcedir data/heart/ \
--prefix heart \
--outdir data/heart

# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/heart/heart_mm10.RData \
--annot data/gencode.vM17.annotation.txt \
--prefix heart \
--outdir data/heart

# Within strain ex vs nonex
Rscript R/diffexpr_within_strain.R \
--counts_matrix data/heart/heart_collapsed_counts_matrix.RData \
--meta_file data/heart/heart-meta-data.txt \
--type non_exercised \
--prefix within-strain-ex-vs-nonex \
--outdir results/heart

# Expression qc
Rscript R/expression_qc.R \
--fpkm_matrix data/heart/heart_collapsed_fpkm_matrix.RData \
--outdir results/heart/qc

# t-SNE sample clustering
Rscript R/clustering.R \
--fpkm_matrix data/heart/heart_collapsed_fpkm_matrix.RData \
--meta_file data/heart/heart-meta-data.txt \
--outdir results/heart/qc

# sample correlation
Rscript R/sample_correlation.R \
--counts_matrix data/heart/heart_collapsed_counts_matrix.RData \
--meta_file data/heart/heart-meta-data.txt \
--outdir results/heart/qc

# expression per sample
Rscript R/expression_per_sample.R \
--fpkm_matrix data/heart/heart_collapsed_fpkm_matrix.RData \
--meta_file data/heart/heart-meta-data.txt \
--title "Heart Specific Genes" \
--prefix heart_specific_genes_expr \
--outdir results/heart/qc \
--gene_list data/heart/heart-specific-genes.txt
