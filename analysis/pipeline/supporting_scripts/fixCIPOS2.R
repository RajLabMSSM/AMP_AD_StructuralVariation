#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# This small function adds CIPOS and CIEND to a vcf file
# Inputs: 
# vcf_file is the path for a vcf file
# output_file is the path for the edited vcf file (it will be compressed)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript fixCIPOS.R input_file.vcf output_file.vcf.gz", call.=FALSE)
} 

vcf_file = args[1]
output_file = args[2]

library(vcfR)
library(tidyr)

vcf = read.vcfR(vcf_file)

cipos = extract.info(vcf,"CIPOS")

if (is.na(cipos[1])){
  # Add CIPOS and CIEND
  info = getINFO(vcf)
  svlen = abs(as.numeric(extract.info(vcf,"SVLEN")))
  ci_int = round(svlen*0.1)
  ci95_int = round(svlen*0.01)
  info2 = paste0("CIPOS=-",ci_int,",",ci_int,";CIEND=-",ci_int,",",ci_int,";CIPOS95=-",ci95_int,",",ci95_int,";CIEND95=-",ci95_int,",",ci95_int,";",info)
  #info2 = paste0("CIPOS=0,0;CIEND=0,0;CIPOS95=0,0;CIEND95=0,0;",info)
  meta = vcf@meta
  ci1 = "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">"
  ci2 = "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">"
  ci3 = "##INFO=<ID=CIEND95,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">"
  ci4 = "##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">"
  meta2 = c(meta[1:which(grepl("##INFO",meta))[1]-1], ci1,ci2,ci3,ci4, meta[which(grepl("##INFO",meta))[1]:length(meta)])
  vcf@fix[,"INFO"] = info2
  vcf@meta = meta2
  write.vcf(x = vcf, file = output_file)
}

cipos95 = extract.info(vcf,"CIPOS95")

if (is.na(cipos95[1])){
  # Add CIPOS95 and CIEND95
  info = getINFO(vcf)
  cipos = extract.info(vcf,"CIPOS")
  cipos = tidyr::separate(data = as.data.frame(cipos), col = cipos, into = c("f1","f2"), sep = ",")
  cipos$f1 = as.numeric(cipos$f1)
  cipos$f2 = as.numeric(cipos$f2)
  cipos_95 = round(cipos*0.1)
  cipos_95 = cipos_95 %>% unite(f_merged, f1, f2, sep=",")
  cipos_95 = cipos_95$f_merged

  ciend = extract.info(vcf,"CIEND")
  ciend = tidyr::separate(data = as.data.frame(ciend), col = ciend, into = c("f1","f2"), sep = ",")
  ciend$f1 = as.numeric(ciend$f1)
  ciend$f2 = as.numeric(ciend$f2)
  ciend_95 = round(ciend*0.1)
  ciend_95 = ciend_95 %>% unite(f_merged, f1, f2, sep=",")
  ciend_95 = ciend_95$f_merged

  info2 = paste0("CIPOS95=",cipos_95,";CIEND95=",ciend_95,";",info)
  #info2 = paste0("CIPOS95=0,0;CIEND95=0,0;",info)
  meta = vcf@meta
  ci3 = "##INFO=<ID=CIEND95,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">"
  ci4 = "##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">"
  meta2 = c(meta[1:which(grepl("##INFO",meta))[1]-1], ci3,ci4, meta[which(grepl("##INFO",meta))[1]:length(meta)])
  vcf@fix[,"INFO"] = info2
  vcf@meta = meta2
  write.vcf(x = vcf, file = output_file)
}


