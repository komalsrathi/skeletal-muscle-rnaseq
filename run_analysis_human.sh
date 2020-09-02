# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/human/human_hg38.RData \
--annot data/gencode.vM17.annotation.txt \
--prefix human \
--outdir data/human

# Within strain exercised vs nonexercised
Rscript R/diffexpr_within_strain.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--meta_file data/human/human-meta-data.txt \
--type 'non_exercised' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix within-strain-ex-vs-nonex \
--excel TRUE \
--text TRUE \
--outdir results/human/within-strain-ex-vs-nonex

