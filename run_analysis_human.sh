# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/human/human_hg38.RData \
--annot data/gencode.v23.annotation.txt \
--prefix human \
--outdir data/human

# 399_2 vs 370-1
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--meta_file data/human/human-meta-data.txt \
--type '399_2, 370_1' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_1' \
--excel TRUE \
--text TRUE \
--outdir results/human/diffexpr-399_2-vs-370_1

# 399_2 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--meta_file data/human/human-meta-data.txt \
--type '399_2, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human/diffexpr-399_2-vs-370_5