# SV-xQTL 

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

