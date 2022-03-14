# Analysis results

## File organization

Results are divided into multiple subdirectories, one for each type of analysis conducted in the project. 

_Note: Description of each file is available within respective subdirectories' README._

The structure of this directory is as follows:  

| Subdirectory | Description |
| :--- | :--- |
| [`SV_calls/`](https://github.com/RajLabMSSM/AMP_AD_StructuralVariation/tree/main/data/SV_calls/) | SV site-frequency data from each study cohort |
| [`LD_tagging/`](https://github.com/RajLabMSSM/AMP_AD_StructuralVariation/tree/main/data/LD_tagging/) | LD between SVs and SNVs |
| [`SV_xQTL/`](https://github.com/RajLabMSSM/AMP_AD_StructuralVariation/tree/main/data/SV_xQTL/) | Full nominal and permuted SV-xQTL summary statistics |
| [`meta_analysis/`](https://github.com/RajLabMSSM/AMP_AD_StructuralVariation/tree/main/data/meta_analysis/) | Full mashR summary statistics for SV-eQTLs across brain regions |
| [`disease_associations/`](https://github.com/RajLabMSSM/AMP_AD_StructuralVariation/tree/main/data/disease_associations/) | Full case-control SV associations summary statistics |

*Individual-level calls and PacBio long-read data can accessed at [Synapse:syn26952206](https://www.synapse.org/#!Synapse:syn26952206)*

## File description

### LD_tagging

Complete LD tagging information between SVs and SNVs mapped in ROSMAP.

| Column | Description |
| :--- | :--- |
| _sv_id_ |  SV ID (ROSMAP) |
| _sv_chr_ | Chromosome of the SV  |
| _sv_pos_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_len_ | SV length |
| _snv_id_ | rsID |
| _snv_chr_ | SNV chromosome |
| _snv_pos_ |  SNV position in the chromosome |
| _r2_ |  LD R2  |

### meta_analysis

Summary statistics from [mashR](https://doi.org/10.1038/s41588-018-0268-8). Meta-analysis was performed across all AMP-AD cohorts and brain regions. Only the best variant per gene is reported. 

| Column | Description |
| :--- | :--- |
| _cohort_region_ |  Cohort and brain region |
| _sv_ |  Merged cohort SV ID |
| _sv_chr_ | Chromosome of the SV  |
| _sv_start_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_end_ | Most downstream genomic position within the chromosome of SV breakpoints |
| _ensembl_ | ENSEMBL_ID of the tested phenotype (gencode v24) |
| _gene_name_ | Gene name |
| _lfsr_ |  Mash local false sign rates |
| _posterior_mean_ |  Mash posterior means  |
| _standard_error_ |  Standard error |
| _mash_beta_ |  Mash slope (posterior mean * standard error)  |

### disease_associations

SV calls from ROS/MAP, Mayo Clinic and MSBB were merged into a combined call set using SURVIVOR68 while requiring 1000 bp maximum distance between breakpoints to merge SVs of the same type. A total of 22,007 SVs identified and all three study groups and with MAF >= 0.01 were selected for the association test. 

Alzheimer's disease status was harmonized across cohorts as previously described:

**ROSMAP**

* AD: cogdx = 4, braaksc >= 4 and ceradsc <= 2
* CONTROL: cogdx = 1, braaksc <= 3 and ceradsc >= 3

**MSBB**

* AD: CDR >= 1, bbscore >= 4 and NP.1 >= 2
* CONTROL: CDR <= 0.5, bbscore <= 3 and NP.1 <= 1

**Mayo Clinic**

* AD: Braak score >= 4 and CERAD > 1
* CONTROL: Braak score <= 3 CERAD < 2

Progressive supranuclear palsy (PSP) case status was obtained from Mayo Clinic pathological diagnosis. 

| Column | Description |
| :--- | :--- |
| _id_ |  Merged cohort SV ID |
| _rosmap_ |  SV ID from ROS/MAP calls  |
| _mayo_ | SV ID from Mayo calls |
| _msbb_ | SV ID from MSBB calls |
| _chrom_ | Chromosome of the SV  |
| _pos_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_end_ | Most downstream genomic position within the chromosome of SV breakpoints |
| _svlen_ | SV size |
| _alt_ | ALT field from VCF |
| _svtype_ |  INFO/SVTYPE from VCF |
| _ac_ |  Allele count  |
| _an_ |  Total number of alleles in called genotypes |
| _ac_hom_ |  Number of homozygous alleles  |
| _ac_het_ |  Number of heterozygous allleles  |
| _maf_ | Minor allele frequency in the merged cohorts |
| _hwe_ | HWE test |
| _description_ | Genomic annotation referring nearest gene (from bumphunter R package) |
| _gene_name_ | Closest gene |
| _slope_ | Slope of association |
| _std_error_ | Standard error |
| _nom_pval_ |  Nominal P-value |
| _fdr_ | FDR adjusted P-value |
| _bonf_ |  Bonferroni adjusted P-value |

### SV_calls

VCF format files for structural variant (SV) calls discovered separately in 1,106 ROS/MAP, 349 Mayo clinic, and 305 MSBB genomes, and comprehending 72,348, 46,197 and 52,451 SVs respectively. Genomic coordinates are aligned against the GRCh37 reference.

### SV-xQTL 

Nominal and permuted summary statistics for each of the molecular phenotype tested: 

* SV-eQTL (gene-mRNA)
* SV-sQTL (splicing junction)
* SV-haQTL (Histone 3 Lysine 9 acetylation, H3K9Ac peaks)
* SV-pQTL (gene-protein)

RNA-seq data for ROS/MAP are from the dorsolateral prefrontal cortex (DLPFC). 

RNA-seq data from MSBB are from four brain regions: BM10 = Brodmann area 10 (part of the frontopolar prefrontal cortex), BM22 = Brodmann area 22 (part of the superior temporal gyrus), BM36 = Brodmann area 36 (part of the fusiform gyrus), and BM44 = Brodmann area 44 (opercular part of the inferior frontal gyrus). 

RNA-seq from Mayo Clinic are from TCX = temporal cortex, CBE = cerebellum.

The ChIP-seq (H3K9Ac) and proteomics data (Tandem mass tag, TMT) are from ROS/MAP cohorts.

Associations were performed using a modified version of [FastQTL](https://github.com/hall-lab/fastqtl), that considers SVs lying within the *cis* window if any part of the spanned region is within the *cis* window **except** inversions, for which one (or both) of the breakpoints must fall are within the cis window. SVs with MAF >= 0.01 were considered for tests. 

###  SV-eQTL 

| Column | Description |
| :--- | :--- |
| _phenotype_ | ID of the tested phenotype (e.g. ENSEMBL_ID for SV-eQTL) |
| _variant_ | ID of the tested SV |
| _lead_variant_ | TRUE if the SV is lead variant for the phenotype after permutation pass |
| _cor_ | Pearson correlation coefficient between genotype dosages and phenotype quantifications |
| _nom_pval_ | The nominal P-value of association |
| _slope_ | The slope associated with the nominal P-value |
| _nvar_ | Number of variants tested in cis for this phenotype |
| _ppval_ | A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype |
| _bpval_ | A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis |
| _sv_chr_ | Chromosome of the SV |
| _sv_start_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_end_ | Most downstream genomic position within the chromosome of SV breakpoints |
| _sv_midpoint_ | Midpoint genomic position within the chromosome between SV breakpoints |
| _svtype_ | SV class |
| _gene_chr_ | Gene chromosome |
| _gene_start_ | Gene genomic start position within the chromosome |
| _gene_end_ | Gene genomic end position within the chromosome |
| _gene_midpoint_ | Gene genomic midpoint position within the chromosome |
| _gene_name_ | Gene symbol |
| _gene_type_ | Gene biotype |
| _dist_midpoints_ | Distance between sv_midpoint and gene_midpoint |
| _dist_closest_breakpoint_ | Distance of the closest SV breakpoint from gene body (if any breakpoint falls within gene body coordinates, distance is equal zero) |
| _cohort_ | Study cohort |
| _total_exon_overlap_bp_ | Total number of phenotype exonic region overlapped by the SV |
| _perc_coding_overlap_ | Percent of phenotype coding region overlapped by the SV |

### SV-sQTL

| Column | Description |
| :--- | :--- |
| _phenotype_ | ID of the tested phenotype (e.g. chr:start-end:cluster_ID for SV-sQTL) |
| _variant_ | ID of the tested SV |
| _lead_variant_ | TRUE if the SV is lead variant for the phenotype after permutation pass |
| _cor_ | Pearson correlation coefficient between genotype dosages and phenotype quantifications |
| _nom_pval_ | The nominal P-value of association |
| _slope_ | The slope associated with the nominal P-value |
| _nvar_ | Number of variants tested in cis for this phenotype |
| _ppval_ | A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype |
| _bpval_ | A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis |
| _sv_chr_ | Chromosome of the SV |
| _sv_start_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_end_ | Most downstream genomic position within the chromosome of SV breakpoints |
| _sv_midpoint_ | Midpoint genomic position within the chromosome between SV breakpoints |
| _svtype_ | SV class |
| _pheno_chr_ | Intron junction chromosome |
| _pheno_start_ | Intron junction genomic upstream position within the chromosome |
| _pheno_end_ | Intron junction genomic downstream position within the chromosome |
| _pheno_name_ | Intron junction cluster ID |
| _gene_chr_ | Gene chromosome |
| _gene_start_ | Gene genomic start position within the chromosome |
| _gene_end_ | Gene genomic end position within the chromosome |
| _gene_midpoint_ | Gene genomic midpoint position within the chromosome |
| _gene_name_ | Gene symbol |
| _dist_midpoints_ | Distance between sv_midpoint and gene_midpoint |
| _dist_closest_breakpoint_ | Distance of the closest SV breakpoint from gene body (if any breakpoint falls within gene body coordinates, distance is equal zero) |
| _dist_closest_breakpoint_pheno_ | Distance of the closest SV breakpoint from the phenotype (i.e. intron junction) |
| _cohort_ | Study cohort |

### SV-haQTL

| Column | Description |
| :--- | :--- |
| _phenotype_ | ID of the tested phenotype (e.g. gene_symbol\|peak_ID\|chr:start-end for SV-haQTL) |
| _variant_ | ID of the tested SV |
| _lead_variant_ | TRUE if the SV is lead variant for the phenotype after permutation pass |
| _cor_ | Pearson correlation coefficient between genotype dosages and phenotype quantifications |
| _nom_pval_ | The nominal P-value of association |
| _slope_ | The slope associated with the nominal P-value |
| _nvar_ | Number of variants tested in cis for this phenotype |
| _ppval_ | A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype |
| _bpval_ | A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis |
| _sv_chr_ | Chromosome of the SV |
| _sv_start_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_end_ | Most downstream genomic position within the chromosome of SV breakpoints |
| _sv_midpoint_ | Midpoint genomic position within the chromosome between SV breakpoints |
| _svtype_ | SV class |
| _pheno_chr_ | Chromosome of the H3K9ac peak |
| _pheno_start_ | H3K9ac peak genomic upstream position within the chromosome |
| _pheno_end_ | H3K9ac peak genomic downstream position within the chromosome |
| _pheno_name_ | H3K9ac peak ID |
| _gene_chr_ | Gene chromosome (closest gene from the peak) |
| _gene_start_ | Gene genomic start position within the chromosome (closest gene from the peak) |
| _gene_end_ | Gene genomic end position within the chromosome (closest gene from the peak) |
| _gene_midpoint_ | Gene genomic midpoint position within the chromosome (closest gene from the peak) |
| _gene_name_ | Gene symbol (closest gene from the peak) |
| _dist_midpoints_ | Distance between sv_midpoint and gene_midpoint |
| _dist_closest_breakpoint_ | Distance of the closest SV breakpoint from gene body (if any breakpoint falls within gene body coordinates, distance is equal zero) |
| _dist_closest_breakpoint_pheno_ | Distance of the closest SV breakpoint from the phenotype (i.e. H3K9ac peak) |
| _cohort_ | Study cohort |

### SV-pQTL

| Column | Description |
| :--- | :--- |
| _phenotype_ | ID of the tested phenotype (e.g. gene_name\|uniprot_id for SV-pQTL) |
| _variant_ | ID of the tested SV |
| _lead_variant_ | TRUE if the SV is lead variant for the phenotype after permutation pass |
| _cor_ | Pearson correlation coefficient between genotype dosages and phenotype quantifications |
| _nom_pval_ | The nominal P-value of association |
| _slope_ | The slope associated with the nominal P-value |
| _nvar_ | Number of variants tested in cis for this phenotype |
| _ppval_ | A first permutation p-value directly obtained from the permutations with the direct method. This is basically a corrected version of the nominal p-value that accounts for the fact that multiple variants are tested per molecular phenotype |
| _bpval_ | A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis |
| _sv_chr_ | Chromosome of the SV |
| _sv_start_ | Most upstream genomic position within the chromosome of SV breakpoints |
| _sv_end_ | Most downstream genomic position within the chromosome of SV breakpoints |
| _sv_midpoint_ | Midpoint genomic position within the chromosome between SV breakpoints |
| _svtype_ | SV class |
| _gene_chr_ | Gene chromosome |
| _gene_start_ | Gene genomic start position within the chromosome |
| _gene_end_ | Gene genomic end position within the chromosome |
| _gene_midpoint_ | Gene genomic midpoint position within the chromosome |
| _gene_name_ | Gene symbol |
| _gene_type_ | Gene biotype |
| _dist_midpoints_ | Distance between sv_midpoint and gene_midpoint |
| _dist_closest_breakpoint_ | Distance of the closest SV breakpoint from gene body (if any breakpoint falls within gene body coordinates, distance is equal zero) |
| _cohort_ | Study cohort |
| _total_exon_overlap_bp_ | Total number of phenotype exonic region overlapped by the SV |
| _perc_coding_overlap_ | Percent of phenotype coding region overlapped by the SV |

