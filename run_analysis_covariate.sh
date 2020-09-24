# mouse skm
echo "Mouse SKM"
Rscript R/gene_cov_analysis.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--covariate_file data/mouse_skm/skm_covariable_list.txt \
--var_filter TRUE \
--col 'label' \
--outdir results/skm/covars

# mouse skm corrected
echo "Mouse SKM corrected"
Rscript R/gene_cov_analysis.R \
--counts_matrix data/mouse_skm/skm_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_skm/skm-meta-data.txt \
--covariate_file data/mouse_skm/skm_covariable_list.txt \
--var_filter TRUE \
--col 'label' \
--outdir results/skm_bc/covars

# mouse heart
echo "Mouse Heart"
Rscript R/gene_cov_analysis.R \
--counts_matrix data/mouse_heart/heart_collapsed_counts_matrix.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--covariate_file data/mouse_heart/heart_covariable_list.txt \
--var_filter TRUE \
--col 'label' \
--outdir results/heart/covars

# mouse heart corrected
echo "Mouse Heart corrected"
Rscript R/gene_cov_analysis.R \
--counts_matrix data/mouse_heart/heart_collapsed_counts_matrix_batchcorrected.RData \
--meta_file data/mouse_heart/heart-meta-data.txt \
--covariate_file data/mouse_heart/heart_covariable_list.txt \
--var_filter TRUE \
--col 'label' \
--outdir results/heart_bc/covars