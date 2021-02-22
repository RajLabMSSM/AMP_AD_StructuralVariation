# Summary statistics results for disease associations

SV calls from ROS/MAP, Mayo Clinic and MSBB were merged into a combined call set using SURVIVOR while requiring 1000bp maximum distance between breakpoints to merge SVs of the same type. 

A total of 22,007 SVs identified and all three study groups and with MAF >= 0.01 were selected for the association test. 

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


