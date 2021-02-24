#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# This small function receives a SURVIVOR merged tools VCF and returns just one genotype for each SV following the rule:
# 1st - Lumpy (svtyper)
# 2nd - Manta
# 3rd - Delly
# 4th - BreakSeq
# 5th - CNVnator
# 6th - BreakDancer (not genotyped) = ./.
# The output vcf will contain only one column with the GT field only. 
#
# Inputs: 
# vcf_file is the path for a vcf file
# sample_name is the id of the sample to be written on the column name
# output_file is the path for the edited vcf file (it will be compressed)

# NOTE: results from this script are not final. SVs are force genotyped later in the pipeline using smoove.

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Usage: Rscript fix_SURVIVORgenotypes.R input_file.vcf sample_name output_file.vcf.gz", call.=FALSE)
}

vcf_file = args[1]
sample_name = args[2]
output_file = args[3]

library(vcfR)

vcf = read.vcfR(vcf_file)

# Get only GT field (must be the 1st in the FORMAT)
#gts = gsub(x = vcf@gt[,-1], pattern = "^(.*?):(.*)", replacement = "\\1")
gts = as.data.frame(matrix(gsub(x = vcf@gt[,-1], pattern = "^(.*?):(.*)", replacement = "\\1"),nrow = nrow(vcf@gt)))
colnames(gts) = colnames(as.data.frame(vcf@gt))[-1]

# Get the column of each tool
manta_idx = grepl(x = colnames(gts),pattern = "MANTA", ignore.case = T)
lumpy_idx = grepl(x = colnames(gts),pattern = "LUMPY", ignore.case = T)
delly_idx = grepl(x = colnames(gts),pattern = "DELLY", ignore.case = T)
breakseq_idx = grepl(x = colnames(gts),pattern = "BREAKSEQ", ignore.case = T)
breakdancer_idx = grepl(x = colnames(gts),pattern = "BREAKDANCER", ignore.case = T)
cnvnator_idx = grepl(x = colnames(gts),pattern = "CNVNATOR", ignore.case = T)

# Rules for genotyping:
# 1st - Lumpy (svtyper)
# 2nd - Manta
# 3rd - Delly
# 4th - BreakSeq
# 5th - CNVnator
# 6th - BreakDancer (not genotyped) = ./.

newgt = matrix("./.", nrow = nrow(gts), ncol = 1)

for (i in 1:nrow(gts)){
  if( sum(gts[i,]!="./." & lumpy_idx)>0 ){
    newgt[i] = as.character(gts[i,lumpy_idx])
    next
  }
  if( sum(gts[i,]!="./." & manta_idx)>0 ){
    newgt[i] = as.character(gts[i,manta_idx])
    next
  }
  if( sum(gts[i,]!="./." & delly_idx)>0 ){
    newgt[i] = as.character(gts[i,delly_idx])
    next
  }
  if( sum(gts[i,]!="./." & breakseq_idx)>0 ){
    newgt[i] = as.character(gts[i,breakseq_idx])
    next
  }
  if( sum(gts[i,]!="./." & cnvnator_idx)>0 ){
    newgt[i] = as.character(gts[i,cnvnator_idx])
    next
  }
  if( sum(gts[i,]!="./." & breakdancer_idx)>0 ){
    newgt[i] = as.character(gts[i,breakdancer_idx])
    next
  }
}

# New FORMAT field will contain only GT
fix_gt = cbind("GT",newgt)
colnames(fix_gt) = c("FORMAT",sample_name)
vcf@gt = fix_gt

# Renaming IDs
info = getINFO(vcf)
sv_type = gsub(pattern = '(.*)(SVTYPE=)(.{3})(;)(.*)', replacement = "\\3", x = info)

sv_type_counts = matrix(seq(0,0, length.out = length(unique(sv_type))))
rownames(sv_type_counts) = unique(sv_type)

for (i in 1:length(sv_type)){
  sv_type_counts[ which(rownames(sv_type_counts) == sv_type[i]) ] = 1 + sv_type_counts[ which(rownames(sv_type_counts) == sv_type[i]) ]
  vcf@fix[i,"ID"] = paste0(sv_type[i],"_",sv_type_counts[ which(rownames(sv_type_counts) == sv_type[i]) ])
}

write.vcf(x = vcf, file = output_file)
