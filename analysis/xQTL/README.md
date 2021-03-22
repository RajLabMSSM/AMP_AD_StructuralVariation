SV-xQTL mapping
================

All tests were performed using a modified version of [FastQTL](https://github.com/hall-lab/fastqtl) to account for SVs. Permuted results were used to determine the lead SV per gene and *P*-values were adjusted for multiple testing using Benjamini-Hochberg (FDR). Associations were performed separately for each SV class. 

Documentation for FastQTL can be found here: http://fastqtl.sourceforge.net

Once the data is processed accordingly. The script [runFastQTL.sh](https://github.com/RajLabMSSM/AMP_AD_StructuralVariation/tree/main/analysis/xQTL/runFastQTL.sh) can be used as an archetype for running the associations for any phenotype.

Details of each phenotype is described below.

------

## SV-eQTL

For SV-eQTL analysis, we used previously processed residuals values obtained after adjusting for clinical, technical, and hidden (via SVA) confounders following conditional quantile normalization (CQN) normalization to account for variations in gene length and GC content and outlier detection and removal. We mapped SV-eQTL for SVs with MAF ≥ 0.01 and within a 1 Mb window from each gene TSS. 

## SV-haQTL

For SV-haQTL analysis, we used residualized values obtained from 571 samples with WGS after regressing out “Sex”, “gel_batch”, “AgeAtDeath” and the first 3 principal components of the genotype matrix to account for the effect of ancestry plus the first 10 principal components of the phenotype matrix to account for the effect of known and hidden factors. We mapped SV-haQTL for SVs with MAF ≥ 0.01 and within a 1 Mb window from each peak. 

## SV-sQTL

For SV-sQTL analysis, we regress out the effects of known and hidden factors as performed in the original [manuscript](https://www.nature.com/articles/s41588-018-0238-1). We accounted for the effect of ancestry given by the first three principal components of the genotype matrix plus the first 15 principal components of the phenotype matrix (PSI). A total of 505 samples with WGS data were used in the SV-sQTL analysis. SVs with MAF ≥ 0.01 and within 100 kb of each intron junction were tested.

## SV-pQTL

For SV-pQTL analysis, we used residualized values for 7,960 protein obtained from 272 samples with WGS after regressing “PMI”, “Sex”, “AgeAtDeath”, three first ancestry PCs, and the first 10 principal components of the phenotype matrix. We tested SVs with MAF ≥ 0.01 and within 1 Mb of each protein.

