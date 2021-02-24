#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

# This small function removes sample outliers in a population VCF file
# Inputs: 
# vcf_file is the path for a vcf file
# output_file is the path for the edited vcf file (it will be compressed)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Usage: Rscript removeSampleOutliers.R input_file.vcf output_file.vcf.gz outlier_samples.blacklist", call.=FALSE)
} 

vcf_file = args[1]
output_file = args[2]
blacklist_file = args[3]

library(vcfR)
library(splitstackshape)
library(abind)
library(reshape2)

######################################################
# Subfunctions from gnomAD-SV pipeline

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.
# (https://github.com/talkowski-lab/gnomad-sv-pipeline/blob/master/gnomad_sv_analysis_scripts/determine_svcount_outliers.R)

#Get list of outliers for a single svtype 
getOutliers <- function(dat,svtype,n.iqr,minPerSVTYPE){
  #Get counts corresponding to svtype
  vals <- dat$count[which(dat$svtype==svtype)]
  if (length(vals)==0)
    return(NA)
  
  #Determine cutoffs
  quartiles <- as.numeric(quantile(vals,probs=c(0.25,0.75),na.rm=T))
  spacer <- n.iqr*IQR(vals,na.rm=T)
  if(median(vals,na.rm=T) >= minPerSVTYPE){
    cutoffs <- c(quartiles[1]-spacer,quartiles[2]+spacer)
  }else{
    cutoffs <- c(NA,NA)
  }
  
  #Return list of sample IDs failing cutoffs
  if(any(!is.na(cutoffs))){
    fails <- dat$sample[which(dat$svtype==svtype
                              & (dat$count<cutoffs[1] | dat$count>cutoffs[2]))]
  }else{
    fails <- as.character(c())
  }
  return(fails)
}

######################################################

get_vcfINFO <- function(vcf){
  info_tbl = read.table(text = (getINFO(vcf)), sep = ';')
  info_tbl2 = as.data.frame(cSplit(info_tbl, 1:ncol(info_tbl), sep = "=", type.convert = FALSE))
  info_tbl3 = data.frame(info_tbl2[,seq(2,ncol(info_tbl2),2)])
  colnames(info_tbl3) = info_tbl2[1,seq(1,ncol(info_tbl2),2)]
  return(info_tbl3)
}

#Defaults
n.iqr = 3
minPerSVTYPE = 50

vcf = read.vcfR(vcf_file)

vcf_sample_names = colnames(as.data.frame(vcf@gt)[,-1])

vcf_info = get_vcfINFO(vcf)

# Select SV type info
vcf_info$SVTYPE = factor(vcf_info$SVTYPE, levels = c("DEL","DUP","INV","INS","TRA","BND"))

gt_only = gsub(pattern = "(.*?):(.*)", replacement = "\\1", x = vcf@gt)
idx_genotyped = gt_only[, 2:ncol(gt_only)] == c("0/1") | gt_only[, 2:ncol(gt_only)] == c("1/1")
idx_genotyped[is.na(idx_genotyped)] = FALSE
vcf_supp_matrix = data.matrix(idx_genotyped)
vcf_info$SUPP = rowSums(vcf_supp_matrix)

#Prepare data from vcf_obj
df = aggregate(vcf_supp_matrix, by = list(vcf_info$SVTYPE), FUN=sum)
# first remember the names
n <- df$Group.1
# transpose all but the first column (name)
df.t <- as.data.frame(t(df[,-1]))
colnames(df.t) <- n
df.t$sample = rownames(df.t)
dat = melt(df.t, id.vars = "sample")
colnames(dat) = c("sample","svtype","count")

outliers = list(DEL = list(),
                DUP = list(),
                INS = list(),
                INV = list(),
                TRA = list(),
                BND = list())

for (svtype in c("DEL","DUP","INS","INV","TRA","BND")){
  outliers[[svtype]] = getOutliers(dat,svtype,n.iqr,minPerSVTYPE)
}

SamplesToRemove = unique(na.omit(unlist(outliers, recursive = T)))
print(paste0(length(SamplesToRemove)," outlier samples excluded."))
if (length(SamplesToRemove)>0){
  write.table(x = data.frame(sample = SamplesToRemove, reason = "sample_outlier"), file = blacklist_file, row.names = F, col.names = F, quote = F, sep = '\t')  
}else{
  system(paste0("cat /dev/null > ", blacklist_file))
}

vcf@gt <- vcf@gt[ , !(colnames(vcf@gt)%in%SamplesToRemove) ]

write.vcf(x = vcf, file = output_file)
