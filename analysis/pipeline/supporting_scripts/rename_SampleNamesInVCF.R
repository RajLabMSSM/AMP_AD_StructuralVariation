#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# This small function renames the vcf sample names
# Inputs: 
# vcf_file is the path for a vcf file (may be compressed)
# new_names is a list of new names separed by semicolon
# output_file is the path for the edited vcf file (it will be compressed)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Usage: Rscript rename_SampleNamesInVCF.R input_file.vcf newname1;newname2 output_file.vcf.gz", call.=FALSE)
} 

vcf_file = args[1]
new_names = args[2]
output_file = args[3]

library(vcfR)

vcf = read.vcfR(vcf_file)

newnames = strsplit(new_names,';')[[1]]
if ( length(newnames) == length(colnames(vcf@gt))-1 ) {
  colnames(vcf@gt)[-1] = newnames
} else {
  stop("Wrong number of ids for this vcf")
}

write.vcf(x = vcf, file = output_file)
