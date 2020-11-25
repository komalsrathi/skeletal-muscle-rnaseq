setwd('~/Projects/Skeletal_Muscle_Patrick/')
library(GSEA)

hallmark_out <- GSEA(input.ds = 'results/skeletal_muscle/gsea/gsea_all_strains.gct',
                     input.cls = 'results/skeletal_muscle/gsea/gsea_all_strains.cls',
                     gs.db = 'data/gsea/h.all.v7.2.symbols.gmt', 
                     nperm = 1000, 
                     random.seed = 42, 
                     gs.size.threshold.min = 15, gs.size.threshold.max = 1500,
                     gsea.type = 'GSEA',
                     weighted.score.type = 1) 

c2_cp_out <- GSEA(input.ds = 'results/skeletal_muscle/gsea/gsea_all_strains.gct',
                  input.cls = 'results/skeletal_muscle/gsea/gsea_all_strains.cls',
                  gs.db = 'data/gsea/c2.cp.v7.2.symbols.gmt', 
                  nperm = 1000, 
                  random.seed = 42, 
                  gs.size.threshold.min = 15, gs.size.threshold.max = 1500,
                  gsea.type = 'GSEA',
                  weighted.score.type = 1) 
