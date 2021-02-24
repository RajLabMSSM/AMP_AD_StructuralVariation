#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 04_MergedSamples
## Usage: ./collect_04_MergedSamples.sh myConfigFile.config

## This script must be ran after all jobs from run_04_MergedSamples.sh are finished.
## A folder Results will be created inside the folder 04_MergedSamples.

#################################################################################

display_usage() {
        echo "Script must run with the configuration file."
        echo -e "\nUsage:\n $0 shortReadPipeline.config \n"
        }

# if less than one argument supplied, display usage
if [[ $# -eq 0 ]] ; then
    display_usage
    exit 1
fi
# check whether user had supplied -h or --help . If yes display usage
if [[ ( $# == "--help") ||  $# == "-h" ]] ; then
    display_usage
    exit 0
fi
# load configFile
source $1
if ! [ -n "$bamFileList" ]; then
    echo "Something is wrong. bamFileList var is unset. Please check if the config file is OK.";
    exit 1
else
    echo "Loading config file seems OK. bamFileList='$bamFileList'.";
    lastline=$(tail -n1 $1)
    if [ "$lastline" = "#### End of file ####" ]; then
        echo "Last line seems OK.";
    else
        echo "Last line must be '#### End of file ####'. Please check if the config file is OK.";
        exit 1
    fi
fi

# =============================================================================
#                           Set folders and modules
# =============================================================================
# Merge all samples (creating a population SV call set of merged tools) using SURVIVOR
##

## Run folder
runDir=${outputFolderPath}"/04_MergedSamples"
mkdir -p ${runDir} || exit 1

ml purge
ml vcftools
ml bcftools
ml samtools
ml bedtools
ml R/3.6.2

callsDir=${outputFolderPath}/03_MergeTools/merged

# =============================================================================
#                              Prepare data
# =============================================================================

# SURVIVOR merging parameters
breakpoint_dist=1000 # max distance between breakpoints
min_num_calls=1 # Minimum number of supporting sample
use_type=1 # Take the type into account (1==yes, else no)
use_strand=1 # Take the strands of SVs into account (1==yes, else no)
dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
min_sv_size=0 # Minimum size of SVs to be taken into account.

# First we select only samples that passed previous filters (i.e. not in any blacklist file)
cat ${outputFolderPath}/sex_check.blacklist \
  ${outputFolderPath}/adapter.blacklist \
  ${outputFolderPath}/aneuploidy.blacklist \
  ${outputFolderPath}/chimera.blacklist \
  ${outputFolderPath}/missing_metadata.blacklist \
  ${outputFolderPath}/svCalls.blacklist | cut -f1 | sort | uniq > ${runDir}/samplesToRemove.list

# Merging samples using merged calls from SURVIVOR 
cat ${id_mapping} | cut -f1 | grep -v -x -f ${runDir}/samplesToRemove.list > ${runDir}/samplesToMerge.list
cat ${id_mapping} | cut -f1 > ${runDir}/samplesToMerge.list

cat /dev/null > ${runDir}/ALL_merged_FileList.list
cat /dev/null > ${runDir}/DEL_merged_FileList.list
cat /dev/null > ${runDir}/DUP_merged_FileList.list
cat /dev/null > ${runDir}/INS_merged_FileList.list
cat /dev/null > ${runDir}/INV_merged_FileList.list
cat /dev/null > ${runDir}/TRA_BND_merged_FileList.list

# Check if sample files exist
while read sample; do
   echo ${callsDir}/${sample}.ALL_merged.vcf >> ${runDir}/ALL_merged_FileList.list
   echo ${callsDir}/${sample}.DEL_merged.vcf >> ${runDir}/DEL_merged_FileList.list
   echo ${callsDir}/${sample}.DUP_merged.vcf >> ${runDir}/DUP_merged_FileList.list
   echo ${callsDir}/${sample}.INS_merged.vcf >> ${runDir}/INS_merged_FileList.list
   echo ${callsDir}/${sample}.INV_merged.vcf >> ${runDir}/INV_merged_FileList.list
   echo ${callsDir}/${sample}.TRA_BND_merged.vcf >> ${runDir}/TRA_BND_merged_FileList.list
   if [ ! -s ${callsDir}/${sample}.ALL_merged.vcf ]; then
      echo "Sample ${sample} failed. Please check"
   fi
done <${runDir}/samplesToMerge.list

# =============================================================================
#                           Merge population call set
# =============================================================================
### Two approaches from here:

# 1 - Merging by SVTYPE first
echo "SURVIVOR merge ${runDir}/DEL_merged_FileList.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/samples_merged_DEL.tmp.vcf" | bsub -n 1 -R 'span[hosts=1]' -R 'rusage[mem=4000]' -P acc_ad-omics -W 12:00 -oo /dev/null -eo /dev/null -q premium
echo "SURVIVOR merge ${runDir}/DUP_merged_FileList.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/samples_merged_DUP.tmp.vcf" | bsub -n 1 -R 'span[hosts=1]' -R 'rusage[mem=4000]' -P acc_ad-omics -W 12:00 -oo /dev/null -eo /dev/null -q premium
echo "SURVIVOR merge ${runDir}/INS_merged_FileList.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/samples_merged_INS.tmp.vcf" | bsub -n 1 -R 'span[hosts=1]' -R 'rusage[mem=4000]' -P acc_ad-omics -W 12:00 -oo /dev/null -eo /dev/null -q premium
echo "SURVIVOR merge ${runDir}/INV_merged_FileList.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/samples_merged_INV.tmp.vcf" | bsub -n 1 -R 'span[hosts=1]' -R 'rusage[mem=4000]' -P acc_ad-omics -W 12:00 -oo /dev/null -eo /dev/null -q premium
echo "SURVIVOR merge ${runDir}/TRA_BND_merged_FileList.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/samples_merged_TRA_BND.tmp.vcf" | bsub -n 1 -R 'span[hosts=1]' -R 'rusage[mem=4000]' -P acc_ad-omics -W 12:00 -oo /dev/null -eo /dev/null -q premium

cat ${runDir}/samples_merged_DEL.tmp.vcf > ${runDir}/samples_merged_ALL.tmp.vcf
grep -v "^#" ${runDir}/samples_merged_DUP.tmp.vcf >> ${runDir}/samples_merged_ALL.tmp.vcf
grep -v "^#" ${runDir}/samples_merged_INS.tmp.vcf >> ${runDir}/samples_merged_ALL.tmp.vcf
grep -v "^#" ${runDir}/samples_merged_INV.tmp.vcf >> ${runDir}/samples_merged_ALL.tmp.vcf
grep -v "^#" ${runDir}/samples_merged_TRA_BND.tmp.vcf >> ${runDir}/samples_merged_ALL.tmp.vcf

bcftools sort -Ov ${runDir}/samples_merged_ALL.tmp.vcf > ${runDir}/samples_merged_ALL.sorted.vcf

${svpipeline_r_lib}/rename_SVIds.R ${runDir}/samples_merged_ALL.sorted.vcf ${runDir}/samples_merged_ALL.sorted.vcf.gz
gunzip -f ${runDir}/samples_merged_ALL.sorted.vcf.gz

# Splitting VCF into BND, DEL, DUP, INS, INV
grep -E '^#|SVTYPE=DEL' ${runDir}/samples_merged_ALL.sorted.vcf > ${runDir}/samples_merged_DEL.sorted.vcf
grep -E '^#|SVTYPE=DUP' ${runDir}/samples_merged_ALL.sorted.vcf > ${runDir}/samples_merged_DUP.sorted.vcf
grep -E '^#|SVTYPE=DEL|SVTYPE=DUP' ${runDir}/samples_merged_ALL.sorted.vcf > ${runDir}/samples_merged_CNV.sorted.vcf
grep -E '^#|SVTYPE=INS' ${runDir}/samples_merged_ALL.sorted.vcf > ${runDir}/samples_merged_INS.sorted.vcf
grep -E '^#|SVTYPE=INV' ${runDir}/samples_merged_ALL.sorted.vcf > ${runDir}/samples_merged_INV.sorted.vcf
grep -E '^#|SVTYPE=BND' ${runDir}/samples_merged_ALL.sorted.vcf > ${runDir}/samples_merged_BND.sorted.vcf
# Delly call TRA as BND. Renaming BND to TRA
sed 's/SVTYPE=BND/SVTYPE=TRA/' ${runDir}/samples_merged_ALL.sorted.vcf | grep -E '^#|SVTYPE=TRA' > ${runDir}/samples_merged_TRA_BND.sorted.vcf

# 2 - Merging all types together and then splitting
SURVIVOR merge ${runDir}/ALL_merged_FileList.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/samples_merged_ALL.vcf

# Rename final SV IDs
${svpipeline_r_lib}/rename_SVIds.R ${runDir}/samples_merged_ALL.vcf ${runDir}/samples_merged_ALL.vcf.gz
gunzip -f ${runDir}/samples_merged_ALL.vcf.gz
cat ${runDir}/samples_merged_ALL.vcf | vcf-sort > ${runDir}/samples_merged_ALL.sorted.vcf
mv ${runDir}/samples_merged_ALL.sorted.vcf ${runDir}/samples_merged_ALL.vcf

# Splitting VCF into BND, DEL, DUP, INS, INV
grep -E '^#|SVTYPE=DEL' ${runDir}/samples_merged_ALL.vcf > ${runDir}/samples_merged_DEL.vcf
grep -E '^#|SVTYPE=DUP' ${runDir}/samples_merged_ALL.vcf > ${runDir}/samples_merged_DUP.vcf
grep -E '^#|SVTYPE=DEL|SVTYPE=DUP' ${runDir}/samples_merged_ALL.vcf > ${runDir}/samples_merged_CNV.vcf
grep -E '^#|SVTYPE=INS' ${runDir}/samples_merged_ALL.vcf > ${runDir}/samples_merged_INS.vcf
grep -E '^#|SVTYPE=INV' ${runDir}/samples_merged_ALL.vcf > ${runDir}/samples_merged_INV.vcf
grep -E '^#|SVTYPE=BND' ${runDir}/samples_merged_ALL.vcf > ${runDir}/samples_merged_BND.vcf
# Delly call TRA as BND. Renaming BND to TRA
sed 's/SVTYPE=BND/SVTYPE=TRA/' ${runDir}/samples_merged_ALL.vcf | grep -E '^#|SVTYPE=TRA' > ${runDir}/samples_merged_TRA_BND.vcf
