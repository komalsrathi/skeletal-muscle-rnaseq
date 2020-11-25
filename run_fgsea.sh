# descending p-value and model r-squared
# skeletal muscle
Rscript --vanilla R/run_fgsea.R \
--input_dir 'results/skeletal_muscle/covars' \
--covars 'delta.running.time, delta.VO2max, RER.avg.CLAMS, respiration.oxphos, running.time, VO2max' \
--type 'd' \
--ranking 'p' \
--output_dir 'results/skeletal_muscle/gsea_preranked'

# heart
Rscript --vanilla R/run_fgsea.R \
--input_dir 'results/heart/covars' \
--covars 'Diameter_d, Stroke.Volume, Ejection.Fraction, Fractional.Shortening, LV.Mass, LVAW_d, LVPW_d' \
--type 'd' \
--ranking 'p' \
--output_dir 'results/heart/gsea_preranked'