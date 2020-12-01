# mouse skeletal muscle
Rscript --vanilla R/create_gsea_input_files.R \
--count_matrix 'data/mouse_skm/skm_collapsed_counts_matrix.RData' \
--meta_file 'data/mouse_skm/skm-meta-data.txt' \
--output_dir 'results/skeletal_muscle/gsea' \
--prefix 'skm'

# mouse heart
Rscript --vanilla R/create_gsea_input_files.R \
--count_matrix 'data/mouse_heart/heart_collapsed_counts_matrix.RData' \
--meta_file 'data/mouse_heart/heart-meta-data.txt' \
--output_dir 'results/heart/gsea' \
--prefix 'heart'
