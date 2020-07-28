# Within strain ex vs nonex
Rscript R/diffexpr_within_strain.R \
--counts_matrix data/collapsed_counts_matrix.RData \
--meta_file data/meta-data-withlitter.txt \
--type non_exercised \
--prefix within-strain-ex-vs-nonex \
--outdir results/