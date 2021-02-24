#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

args = commandArgs(trailingOnly=TRUE)

if (!length(args)%in%c(2,3)) {
  stop("Usage: Rscript complementAnnotation.R anotFile.txt output_file.txt ref", call.=FALSE)
} 

annot_input = args[1]
output_file = args[2]
ref = args[3]

library(bumphunter)
library(dplyr)

if (ref=="h38"){
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
  genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
}else{
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
}

annot = read.delim2(annot_input, sep = '\t', header = T, stringsAsFactors = F)
tab <- matchGenes(makeGRangesFromDataFrame(annot[,c("SV.chrom","SV.start","SV.end")] %>% mutate(SV.chrom = paste0("chr",SV.chrom)) ),genes)
annot = cbind(annot, tab)

write.table(annot, file = output_file, sep = "\t", quote = F, row.names = F, col.names = T)
