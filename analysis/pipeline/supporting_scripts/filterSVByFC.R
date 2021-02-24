#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript filterSVByFC.R input_file.vcf output_file.vcf.gz", call.=FALSE)
} 

vcf_file = args[1]
output_file = args[2]

library(vcfR)
library(tidyr)
library(ggsci)
library(ggplot2)
library(gridExtra)

vcf = read.vcfR(vcf_file, nrows = 1, verbose = F)

# Check if VCF contains the FORMAT fields DHFFC and DHBFC
ffc = grep(pattern = "DHFFC", x = vcf@gt[1])
bfc = grep(pattern = "DHBFC", x = vcf@gt[1])
if (is.na(ffc[1]) | is.na(bfc[1])){
  stop("VCF file needs Duphold FORMAT fields DHFFC and DHBFC", call.=FALSE)
}
# If its OK, proceed
vcf = read.vcfR(vcf_file)

vcf.gt = as.data.frame(extract.gt(vcf,"GT"))
vcf.ffc = as.data.frame(extract.gt(vcf,"DHFFC"))
vcf.bfc = as.data.frame(extract.gt(vcf,"DHBFC"))

indx <- sapply(vcf.ffc, is.factor)
vcf.ffc[indx] <- lapply(vcf.ffc[indx], function(x) as.numeric(as.character(x)))

indx <- sapply(vcf.bfc, is.factor)
vcf.bfc[indx] <- lapply(vcf.ffc[indx], function(x) as.numeric(as.character(x)))

# Read SVTYPES
vcf.svtype = data.frame(SVTYPE = extract.info(vcf,"SVTYPE"))

# Use the thresholds of DHFFC < 0.7 for DEL and DHBFC > 1.3 for DUP as suggested by Brent Pedersen (https://github.com/brentp/duphold)
## DEL
DEL_comp = data.frame( num_ind_with_DEL_in_threshold = rowSums( (vcf.gt=="0/1" | vcf.gt=="1/1") & vcf.ffc<0.7 , na.rm = T),
                       num_ind_with_DEL_gt = rowSums( (vcf.gt=="0/1" | vcf.gt=="1/1") , na.rm = T))
DEL_comp$toKeep = DEL_comp$num_ind_with_DEL_in_threshold/DEL_comp$num_ind_with_DEL_gt>0.7
DEL_to_keep = DEL_comp$toKeep==T & vcf.svtype$SVTYPE=="DEL"

## DUP
DUP_comp = data.frame( num_ind_with_DUP_in_threshold = rowSums( (vcf.gt=="0/1" | vcf.gt=="1/1") & vcf.bfc>1.3 , na.rm = T),
                       num_ind_with_DUP_gt = rowSums( (vcf.gt=="0/1" | vcf.gt=="1/1") , na.rm = T))
DUP_comp$toKeep = DUP_comp$num_ind_with_DUP_in_threshold/DUP_comp$num_ind_with_DUP_gt>0.7
DUP_to_keep = DUP_comp$toKeep==T & vcf.svtype$SVTYPE=="DUP"

png(filename = "samples_filtered_by_FC.png", res = 300)
a = ggplot(DEL_comp[vcf.svtype$SVTYPE=="DEL",], aes(x = num_ind_with_DEL_in_threshold, y = num_ind_with_DEL_gt, color = toKeep)) +
  geom_point(alpha = .3) + theme_classic() + scale_color_d3() + labs(title = "DEL", x="Carrier samples within threshold", y="Carrier samples")

b = ggplot(DUP_comp[vcf.svtype$SVTYPE=="DUP",], aes(x = num_ind_with_DUP_in_threshold, y = num_ind_with_DUP_gt, color = toKeep)) +
  geom_point(alpha = .3) + theme_classic() + scale_color_d3() + labs(title = "DUP", x="Carrier samples within threshold", y="Carrier samples")

grid.arrange(a,b)
dev.off()

vcf_filt = vcf[DEL_to_keep | DUP_to_keep | !vcf.svtype$SVTYPE %in% c("DEL","DUP")]
print(paste0("Variants filtered = ", nrow(vcf@fix) - nrow(vcf_filt@fix)))

write.vcf(x = vcf_filt, file = output_file)
