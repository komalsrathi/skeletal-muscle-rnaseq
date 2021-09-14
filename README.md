# Mouse RNA-Seq Analysis

Authors: Komal S. Rathi (@komalsrathi)

## Introduction

Various steps in the analysis:

1. QC
2. Process using STAR-RSEM pipeline
3. Differential gene expression
4. Pathway analysis using GSEA
5. Regression analysis to identify association between gene expression and clinical covariates
6. Regression analysis to identify association between clinical covariates and strain
7. fGSEA using covariable output
8. Figures

## Methods

### Skeletal Muscle

60 mouse skeletal muscle RNA-sequencing fastq files, consisting of 6 non-exercised, 3 exercised responder and 3 exercised non-responder samples for each of the five strains i.e. ANT1, B6, B6 no ND5, EC77 and IAI, were processed using the STAR alignment tool and subsequently normalized using the RSEM package based upon the mm10 reference genome and the gencode version M17 gene annotation. 

QC was performed by clustering all samples using Principal Component Analysis (PCA). 

Differential gene expression analysis was performed by comparing each gene in the exercised vs the non-exercised group within each strain. The voom procedure was used to normalize the RSEM generated expected counts followed by differential expression testing using R package limma to obtain P-values and LogFC. Specifically, a total of 58581 genes were tested for differential expression between the control and treatment samples. Pathway enrichment was performed using Gene Set Enrichment Analysis (GSEA) version 4.1.0 using a weighted scoring scheme and Hallmark and C2 CP genesets. 

The same procedure of differential gene expression and pathway enrichment using GSEA was repeated for comparisons of 1) non-exercised samples between each strain, 2) exercised samples between each strain and 3) exercised responders vs exercised non-responders across all strains.

Regression analysis was performed in order to find gene expression profiles associated with clinical covariables within each strain and across all strains. Similarly, linear regression was also performed to find covariables associated with strain.

Next, we wanted to get pathways enriched for the genes associated with select covariables like delta running time and delta VO2max for exercised mice in each of the 5 strains as well as across all strains. The covariate analysis output was first split into positively correlated genes with estimate > 0 and negatively correlated genes with estimate < 0. Next, both lists were sorted by ascending p-values and descending r-squared values. Sequential ranks were assigned from most negatively correlated to the most-positively correlated genes from low to high. Using the ranked covariate output, fGSEA (fast preranked gene set enrichment analysis) analysis was performed using the R package fgsea. This process was repeated with covariables like delta running time, delta VO2max, avg CLAMS, respiration oxphos, running time and VO2max for each of the 5 strains as well as across all strains.

### Heart

60 mouse heart RNA-sequencing fastq files, consisting of 6 non-exercised, 6 exercised samples for each of the five strains i.e. ANT1, B6, B6 no ND5, EC77 and IAI, were processed using the STAR alignment tool and subsequently normalized using the RSEM package based upon the mm10 reference genome and the gencode version M17 gene annotation.

QC was performed by clustering all samples using Principal Component Analysis (PCA). 

Differential gene expression analysis was performed by comparing each gene in the exercised vs the non-exercised group within each strain. The voom procedure was used to normalize the RSEM generated expected counts followed by differential expression testing using R package limma to obtain P-values and LogFC. Specifically, a total of 58581 genes were tested for differential expression between the control and treatment samples. Pathway enrichment was performed using Gene Set Enrichment Analysis (GSEA) version 4.1.0 using a weighted scoring scheme and Hallmark and C2 CP genesets. 

The same procedure of differential gene expression and pathway enrichment using GSEA was repeated for comparisons of 1) non-exercised samples between each strain and 2) exercised samples between each strain.

Regression analysis was performed in order to find gene expression profiles associated with clinical covariables within each strain and across all strains. Similarly, linear regression was also performed to find covariables associated with strain.

Next, we wanted to get pathways enriched for the genes associated with select covariables like Diameter, Stroke Volume, Ejection Fraction, Fractional Shortening, LV Mass, LVAW_d and LVPW_d for each of the 5 strains as well as across all strains. The covariate analysis output was first split into positively correlated genes with estimate > 0 and negatively correlated genes with estimate < 0. Next, both lists were sorted by ascending p-values and descending r-squared values. Sequential ranks were assigned from most negatively correlated to the most-positively correlated genes from low to high. Using the ranked covariate output, fGSEA (fast preranked gene set enrichment analysis) analysis was performed using the R package fgsea.

