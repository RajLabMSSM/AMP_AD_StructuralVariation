#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# Inputs: 
# vcf_file is the path for a vcf file
# output_file is the path for the edited vcf file (it will be compressed)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript filterINS.R input_file.vcf output_file.vcf.gz", call.=FALSE)
} 

vcf_file = args[1]
output_file = args[2]

library(vcfR)

vcf = read.vcfR(vcf_file, verbose = F)

fix = as.data.frame(getFIX(vcf))
info = vcfR::extract_info_tidy(vcf)

# We will remove all SVs with SVLEN<50 or ALT == <BND>
svlen = info$SVLEN
alt = fix$ALT

toRemove = svlen<50 | alt=="<BND>"
  
write.vcf(x = vcf[!toRemove,], file = output_file)
