#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript harmonize_insertions cohort_root_dir output_vcf", call.=FALSE)
} 

work_dir = args[1]
output_vcf = args[2]

library(vcfR)
library(GenomicRanges)
library(regioneR)

getoverlap <- function(ref,test){
  refGR <- toGRanges(ref[,1:5])
  testGR <- toGRanges(test[,1:5])
  
  hits = findOverlaps(refGR, testGR)
  overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
  refPercentOverlap <- width(overlaps) / width(refGR[queryHits(hits)])
  testPercentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
  
  hits <- hits[(refPercentOverlap > 0.7) & (testPercentOverlap > 0.7)]
  return(hits)
}

## INS from other callers
survivor_dir = paste0(work_dir, "/03_MergeTools/SURVIVOR/population")
ins = read.vcfR(paste0(survivor_dir,"/samples_merged_INS.vcf"), verbose = F)

ins_bed = as.data.frame(cbind(ins@fix[,c("CHROM","POS","ID")], extract.info(ins, "SVLEN")), stringsAsFactors = F)
colnames(ins_bed) = c("CHROM","POS","ID","SVLEN")
ins_bed$POS = as.numeric(ins_bed$POS)
ins_bed$END = ins_bed$POS + abs(as.numeric(ins_bed$SVLEN))
ins_bed = ins_bed[,c("CHROM","POS","END","ID","SVLEN")]
ins_bed = ins_bed[ order(ins_bed[,1],ins_bed[,2],ins_bed[,3]) , ]

## INS from MELT
melt_dir = paste0(work_dir, "/02_SVtools/MELT/Results")
alu = read.vcfR(paste0(melt_dir,"/ALU.final_comp.vcf"), verbose = F)
line1 = read.vcfR(paste0(melt_dir,"/LINE1.final_comp.vcf"), verbose = F)
sva = read.vcfR(paste0(melt_dir,"/SVA.final_comp.vcf"), verbose = F)

names_alu = colnames(alu@gt)
names_line1 = colnames(line1@gt)
names_sva = colnames(sva@gt)

colnames(alu@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(alu@gt))
colnames(line1@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(line1@gt))
colnames(sva@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(sva@gt))

alu.pass = alu[alu@fix[,7]=="PASS"]
line1.pass = line1[line1@fix[,7]=="PASS"]
sva.pass = sva[sva@fix[,7]=="PASS"]

alu_bed = as.data.frame(cbind(alu.pass@fix[,c("CHROM","POS","ALT")], extract.info(alu.pass, "SVLEN")), stringsAsFactors = F)
colnames(alu_bed) = c("CHROM","POS","ID","SVLEN")
alu_bed$POS = as.numeric(alu_bed$POS)
alu_bed$END = alu_bed$POS + as.numeric(alu_bed$SVLEN)
alu_bed = alu_bed[,c("CHROM","POS","END","ID","SVLEN")]
alu_bed = alu_bed[ order(alu_bed[,1],alu_bed[,2],alu_bed[,3]) , ]

line1_bed = as.data.frame(cbind(line1.pass@fix[,c("CHROM","POS","ALT")], extract.info(line1.pass, "SVLEN")), stringsAsFactors = F)
colnames(line1_bed) = c("CHROM","POS","ID","SVLEN")
line1_bed$POS = as.numeric(line1_bed$POS)
line1_bed$END = line1_bed$POS + as.numeric(line1_bed$SVLEN)
line1_bed = line1_bed[,c("CHROM","POS","END","ID","SVLEN")]
line1_bed = line1_bed[ order(line1_bed[,1],line1_bed[,2],line1_bed[,3]) , ]

sva_bed = as.data.frame(cbind(sva.pass@fix[,c("CHROM","POS","ALT")], extract.info(sva.pass, "SVLEN")), stringsAsFactors = F)
colnames(sva_bed) = c("CHROM","POS","ID","SVLEN")
sva_bed$POS = as.numeric(sva_bed$POS)
sva_bed$END = sva_bed$POS + as.numeric(sva_bed$SVLEN)
sva_bed = sva_bed[,c("CHROM","POS","END","ID","SVLEN")]
sva_bed = sva_bed[ order(sva_bed[,1],sva_bed[,2],sva_bed[,3]) , ]

ins_alu = getoverlap(ins_bed,alu_bed)
ins_line1 = getoverlap(ins_bed,line1_bed)
ins_sva = getoverlap(ins_bed,sva_bed)

# Flag IDs to remove
ins_to_remove = ins_bed[unique(c(queryHits(ins_alu),queryHits(ins_line1),queryHits(ins_sva))),]
write.table(x = ins_to_remove$ID, file = paste0(survivor_dir,"/INS.toremove"), quote = F, row.names = F, col.names = F)

ins_clean = ins[!ins_bed$ID%in%ins_to_remove$ID]
samples = colnames(ins@gt)[-1]

alu.pass@gt = alu.pass@gt[,c("FORMAT",samples)]
line1.pass@gt = line1.pass@gt[,c("FORMAT",samples)]
sva.pass@gt = sva.pass@gt[,c("FORMAT",samples)]

merged = ins_clean

meta = merged@meta
meta2 = c(
  '##fileformat=VCFv4.1',
  paste0(meta[which(grepl("##source",meta))],";",gsub("##source=","",alu.pass@meta[which(grepl("##source",alu.pass@meta))])),
  paste0(meta[which(grepl("##fileDate",meta))],
         ";",gsub("##fileDate=","",alu.pass@meta[which(grepl("##fileDate",alu.pass@meta))]),
         ";",gsub("##fileDate=","",line1.pass@meta[which(grepl("##fileDate",line1.pass@meta))]),
         ";",gsub("##fileDate=","",sva.pass@meta[which(grepl("##fileDate",sva.pass@meta))]))
  )
meta3 = c(meta2,
  c('##reference=human_g1k_v37.fasta',
  '##contig=<ID=1,length=249250621>',
  '##contig=<ID=2,length=243199373>',
  '##contig=<ID=3,length=198022430>',
  '##contig=<ID=4,length=191154276>',
  '##contig=<ID=5,length=180915260>',
  '##contig=<ID=6,length=171115067>',
  '##contig=<ID=7,length=159138663>',
  '##contig=<ID=8,length=146364022>',
  '##contig=<ID=9,length=141213431>',
  '##contig=<ID=10,length=135534747>',
  '##contig=<ID=11,length=135006516>',
  '##contig=<ID=12,length=133851895>',
  '##contig=<ID=13,length=115169878>',
  '##contig=<ID=14,length=107349540>',
  '##contig=<ID=15,length=102531392>',
  '##contig=<ID=16,length=90354753>',
  '##contig=<ID=17,length=81195210>',
  '##contig=<ID=18,length=78077248>',
  '##contig=<ID=19,length=59128983>',
  '##contig=<ID=20,length=63025520>',
  '##contig=<ID=21,length=48129895>',
  '##contig=<ID=22,length=51304566>',
  '##contig=<ID=X,length=155270560>',
  '##contig=<ID=Y,length=59373566>',
  '##ALT=<ID=DEL,Description="Deletion">',
  '##ALT=<ID=DUP,Description="Duplication">',
  '##ALT=<ID=INV,Description="Inversion">',
  '##ALT=<ID=BND,Description="Translocation">',
  '##ALT=<ID=INS,Description="Insertion">',
  '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
  '##ALT=<ID=INS:ME:LINE1,Description="Insertion of LINE1 element">',
  '##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">',
  '##INFO=<ID=CIEND,Number=2,Type=String,Description="PE confidence interval around END">',
  '##INFO=<ID=CIPOS,Number=2,Type=String,Description="PE confidence interval around POS">',
  '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
  '##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">',
  '##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">',
  '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
  '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">',
  '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">',
  '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method for generating this merged VCF file.">',
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">',
  '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of supporting samples.">',
  '##INFO=<ID=SUPP,Number=1,Type=String,Description="Number of samples supporting the variant">',
  '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Indicating the direction of the reads with respect to the type and breakpoint.">',
  '##INFO=<ID=ASSESS,Number=1,Type=Integer,Description="Provides information on evidence availible to decide insertion site.0 = No overlapping reads at site;1 = Imprecise breakpoint due to greater than expected distance between evidence;2 = discordant pair evidence only -- No split read information;3 = left side TSD evidence only;4 = right side TSD evidence only;5 = TSD decided with split reads, highest possible quality">',
  '##INFO=<ID=TSD,Number=1,Type=String,Description="Precise Target Site Duplication for bases, if unknown, value will be NULL">',
  '##INFO=<ID=INTERNAL,Number=2,Type=String,Description="If insertion internal or close to a gene, listed here followed by a discriptor of the location in the gene (either INTRON, EXON_#, 5_UTR, 3_UTR, PROMOTER, or TERMINATOR). If multiple genes intersected, will be seperated by pipe">',
  '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY; If START or END is unknown, will be -1; If POLARITY is unknown, will be null">',
  '##INFO=<ID=DIFF,Number=.,Type=String,Description="Coverage and Differences in relation to the MEI reference. Form is %2XCoverage:Differences, with differences delimited by ,">',
  '##INFO=<ID=LP,Number=1,Type=Integer,Description="Total number of discordant pairs supporting the left side of the breakpont">',
  '##INFO=<ID=RP,Number=1,Type=Integer,Description="Total number of discordant pairs supporting the right side of the breakpont">',
  '##INFO=<ID=RA,Number=1,Type=Float,Description="Ratio between LP and RP, reported as log2(LP / RP)">',
  '##INFO=<ID=PRIOR,Number=1,Type=String,Description="True if this site was not discovered in this dataset, but was included on a provided priors list">',
  '##INFO=<ID=SR,Number=1,Type=Integer,Description="Total number of SRs at the estimated breakpoint for this site. Recomended to filter sites with <= 2 SRs">',
  '##FILTER=<ID=s25,Description="Greater than 25.0% of samples do not have data">',
  '##FILTER=<ID=rSD,Description="Ratio of LP to RP is greater than 2.0 standard deviations">',
  '##FILTER=<ID=hDP,Description="More than the expected number of discordant pairs at this site are also split">',
  '##FILTER=<ID=ac0,Description="No individuals in this VCF file were identified with this insertion">',
  '##FILTER=<ID=lc,Description="MEI is embeded in a low complexity region">',
  '##FILTER=<ID=PASS,Description="Passed all filters">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '##FORMAT=<ID=PSV,Number=1,Type=String,Description="Previous support vector">',
  '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="predicted length">',
  '##FORMAT=<ID=DR,Number=2,Type=Integer,Description="# supporting reference,variant reads in that order">',
  '##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand of SVs">',
  '##FORMAT=<ID=QV,Number=1,Type=String,Description="Quality values: if not defined a . otherwise the reported value.">',
  '##FORMAT=<ID=TY,Number=1,Type=String,Description="Types">',
  '##FORMAT=<ID=ID,Number=1,Type=String,Description="Variant ID from input.">',
  '##FORMAT=<ID=RAL,Number=1,Type=String,Description="Reference allele sequence reported from input.">',
  '##FORMAT=<ID=AAL,Number=1,Type=String,Description="Alternative allele sequence reported from input.">',
  '##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinates">')
)
merged@meta = meta3
merged@gt = rbind(gsub("\\./\\.","0/0",ins_clean@gt), alu.pass@gt, line1.pass@gt, sva.pass@gt)
merged@fix = rbind(merged@fix, alu.pass@fix, line1.pass@fix, sva.pass@fix)
merged@fix[,8] = paste0("SVTYPE=",extract.info(merged,"SVTYPE"),";","SVLEN=",extract.info(merged,"SVLEN"))
merged@gt = gsub("(.*?):(.*)","\\1",merged@gt)

write.vcf(merged, file = output_vcf)
