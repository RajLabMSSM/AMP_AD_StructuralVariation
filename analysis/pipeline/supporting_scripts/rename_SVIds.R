#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# This small function renames variant IDs to $SVTYPE_Number
# Inputs: 
# vcf_file is the path for a vcf file
# output_file is the path for the edited vcf file (it will be compressed)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript rename_SVIds.R input_file.vcf output_file.vcf.gz", call.=FALSE)
} 

vcf_file = args[1]
output_file = args[2]

library(vcfR)

vcf = read.vcfR(vcf_file)

info = getINFO(vcf)
if (sum(grepl("SVTYPE", info))==nrow(as.data.frame(info))){
  sv_type = gsub(pattern = '(.*)(SVTYPE=)(.{3})(.*)', replacement = "\\3", x = info)
}else{
  fix = as.data.frame(getFIX(vcf,getINFO = T),stringsAsFactors = F)
  sv_type = substr(gsub(pattern = '<(.*?)>(.*)', replacement = "\\1", x = fix$ALT), 1, 3)
  vcf@fix[,"INFO"] = paste0("SVTYPE=",sv_type,";",vcf@fix[,"INFO"])
}

sv_type_counts = matrix(seq(0,0, length.out = length(unique(sv_type))))
rownames(sv_type_counts) = unique(sv_type)

for (i in 1:length(sv_type)){
  sv_type_counts[ which(rownames(sv_type_counts) == sv_type[i]) ] = 1 + sv_type_counts[ which(rownames(sv_type_counts) == sv_type[i]) ]
  vcf@fix[i,"ID"] = paste0(sv_type[i],"_",sv_type_counts[ which(rownames(sv_type_counts) == sv_type[i]) ])
}

write.vcf(x = vcf, file = output_file)
