#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# This small function renames the vcf fields REF and ALT to "N" and <SVTYPE>
# Inputs: 
# vcf_file is the path for a vcf file
# output_file is the path for the edited vcf file (it will be compressed)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript rename_REFandALT.R input_file.vcf output_file.vcf.gz", call.=FALSE)
} 

vcf_file = args[1]
output_file = args[2]

library(vcfR)

vcf = read.vcfR(vcf_file)

info = getINFO(vcf)
sv_type = gsub(pattern = '(.*)(SVTYPE=)(.{3})(;)(.*)', replacement = "\\3", x = info)

vcf@fix[,"REF"] = "N"
OLD_ALT = vcf@fix[,"ALT"]
vcf@fix[,"ALT"] = paste0("<",sv_type,">")
vcf@fix[which(sv_type=="INS"),"ALT"] = OLD_ALT[which(sv_type=="INS")]

write.vcf(x = vcf, file = output_file)
