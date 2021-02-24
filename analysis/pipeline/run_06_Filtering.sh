#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 06_Filtering
## Usage: ./run_06_Filtering.sh myConfigFile.config

## A folder named "06_Filtering" will be created

#################################################################################

display_usage() {
        echo "Script must run with the configuration file."
        echo -e "\nUsage:\n $0 myConfigFile.config \n"
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

mkdir -p ${outputFolderPath} || exit 1
runDir=${outputFolderPath}"/06_Filtering"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1
mkdir -p ${runDir}/fromSMOOVE || exit 1
mkdir -p ${runDir}/fromMELT || exit 1

sd=~/ad-omics/ricardo/Data/1000G_phase1/segmental_dups.hg19.bed.gz
sr=~/ad-omics/ricardo/Data/1000G_phase1/simpleRepeats.hg19.bed.gz

## Results from previous modules
genotypesDir_SM=${outputFolderPath}"/05_GenotypeCalls" # (Smoove)

ml purge
ml R
ml vcftools
ml bcftools
ml parallel
ml samtools
ml picard
ml bedtools

# =============================================================================
#                   Harmonizing INS MELT genotyped calls
# =============================================================================

# MELT calls are merged already. So, here we will harmonize those results with other insertion calls.
${svpipeline_r_lib}/harmonize_insertions.R ${outputFolderPath} ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz

#cp ${runDir}/fromSURVIVOR/INS.merged.raw.vcf.gz ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz
${svpipeline_r_lib}/filterINS.R ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz ${runDir}/fromMELT/samples_merged_INS.raw2.vcf.gz
zcat ${runDir}/fromMELT/samples_merged_INS.raw2.vcf.gz | bcftools sort -Oz -o ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz 

# Remove sample outliers
${svpipeline_r_lib}/removeSampleOutliers.R ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz ${runDir}/fromMELT/samples_merged_INS.filt1.vcf.gz ${runDir}/fromMELT/outlier_samples.INS.blacklist

# Remove blacklisted samples again
cat  ${outputFolderPath}/adapter.blacklist \
  ${outputFolderPath}/aneuploidy.blacklist \
  ${outputFolderPath}/chimera.blacklist \
  ${outputFolderPath}/missing_metadata.blacklist \
  ${outputFolderPath}/outlier_samples.blacklist \
  ${outputFolderPath}/sex_check.blacklist \
  ${outputFolderPath}/svCalls.blacklist \
  ${outputFolderPath}/missing_metadata.blacklist \
  ${outputFolderPath}/sex_check.blacklist \
  ${outputFolderPath}/svCalls.blacklist | cut -f1 | sort | uniq > ${runDir}/fromMELT/samplesToRemove.list
bcftools query -l ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz | grep -v -x -f ${runDir}/fromMELT/samplesToRemove.list > ${runDir}/fromMELT/samplesToMerge.list

# Remove blacklisted samples
zcat ${runDir}/fromMELT/samples_merged_INS.raw.vcf.gz | bcftools view --force-samples -S ${runDir}/fromMELT/samplesToMerge.list -Oz > ${runDir}/fromMELT/samples_merged_INS.filt2.vcf.gz

# Filter variants with no ALT allele in any sample. Filter for call rate. Sort variants and clean header.
bcftools view -c1 ${runDir}/fromMELT/samples_merged_INS.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o ${runDir}/fromMELT/samples_merged_INS.filt3.vcf.gz
rm out.log

# Rename IDs
${svpipeline_r_lib}/rename_SVIds.R ${runDir}/fromMELT/samples_merged_INS.filt3.vcf.gz ${runDir}/fromMELT/samples_merged_INS.filt4.vcf.gz

# Copy Final Results
zcat ${runDir}/fromMELT/samples_merged_INS.filt4.vcf.gz | bcftools +fill-tags -Oz -o ${runDir}/fromMELT/samples_merged_INS.Final.vcf.gz

# Filter by MAF
bcftools view -q 0.05:minor ${runDir}/fromMELT/samples_merged_INS.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromMELT/samples_merged_INS.maf05.Final.vcf.gz
bcftools view -q 0.01:minor ${runDir}/fromMELT/samples_merged_INS.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromMELT/samples_merged_INS.maf01.Final.vcf.gz

parallel tabix -f -p vcf {} ::: ${runDir}/fromMELT/*.Final.vcf.gz

# Get tables for plots
bcftools view ${runDir}/fromMELT/samples_merged_INS.Final.vcf.gz | bcftools +missing2ref | bcftools +fill-tags | bcftools sort | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' > ${runDir}/fromMELT/INS.dat
sed -e 's/^/chr/' ${runDir}/fromMELT/INS.dat > ${runDir}/fromMELT/INS.chr.dat
bedtools annotate -both -i ${runDir}/fromMELT/INS.chr.dat -files ${sd} ${sr} > ${runDir}/fromMELT/INS.sd_sr_cov.txt

# =============================================================================
#                   Merging and filtering Smoove genotyped calls
# =============================================================================

# Copy population call after Smoove (merged with Smoove)
cp ${genotypesDir_SM}/Results/DEL.smoove.square.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DEL.raw.vcf.gz
cp ${genotypesDir_SM}/Results/DUP.smoove.square.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DUP.raw.vcf.gz
cp ${genotypesDir_SM}/Results/INV.smoove.square.vcf.gz ${runDir}/fromSMOOVE/samples_merged_INV.raw.vcf.gz

# Remove sample outliers
${svpipeline_r_lib}/removeSampleOutliers.R ${runDir}/fromSMOOVE/samples_merged_DEL.raw.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DEL.filt1.vcf.gz ${runDir}/fromSMOOVE/outlier_samples.DEL.SM.blacklist
${svpipeline_r_lib}/removeSampleOutliers.R ${runDir}/fromSMOOVE/samples_merged_DUP.raw.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DUP.filt1.vcf.gz ${runDir}/fromSMOOVE/outlier_samples.DUP.SM.blacklist
${svpipeline_r_lib}/removeSampleOutliers.R ${runDir}/fromSMOOVE/samples_merged_INV.raw.vcf.gz ${runDir}/fromSMOOVE/samples_merged_INV.filt1.vcf.gz ${runDir}/fromSMOOVE/outlier_samples.INV.SM.blacklist

cat ${runDir}/fromSMOOVE/outlier_samples.DEL.SM.blacklist \
  ${runDir}/fromSMOOVE/outlier_samples.DUP.SM.blacklist \
  ${runDir}/fromSMOOVE/outlier_samples.INV.SM.blacklist | sort | uniq > ${outputFolderPath}/outlier_samples.blacklist

# Remove blacklisted samples again
cat  ${outputFolderPath}/adapter.blacklist \
  ${outputFolderPath}/aneuploidy.blacklist \
  ${outputFolderPath}/chimera.blacklist \
  ${outputFolderPath}/missing_metadata.blacklist \
  ${outputFolderPath}/outlier_samples.blacklist \
  ${outputFolderPath}/sex_check.blacklist \
  ${outputFolderPath}/svCalls.blacklist \
  ${outputFolderPath}/missing_metadata.blacklist \
  ${outputFolderPath}/sex_check.blacklist \
  ${outputFolderPath}/svCalls.blacklist | cut -f1 | sort | uniq > ${runDir}/fromSMOOVE/samplesToRemove.list

cat ${id_mapping} | cut -f1 | grep -v -x -f ${runDir}/fromSMOOVE/samplesToRemove.list > ${runDir}/fromSMOOVE/samplesToMerge.list

# Remove blacklisted samples
bcftools view --force-samples -S ${runDir}/fromSMOOVE/samplesToMerge.list ${runDir}/fromSMOOVE/samples_merged_DEL.raw.vcf.gz -Oz > ${runDir}/fromSMOOVE/samples_merged_DEL.filt2.vcf.gz
bcftools view --force-samples -S ${runDir}/fromSMOOVE/samplesToMerge.list ${runDir}/fromSMOOVE/samples_merged_DUP.raw.vcf.gz -Oz > ${runDir}/fromSMOOVE/samples_merged_DUP.filt2.vcf.gz
bcftools view --force-samples -S ${runDir}/fromSMOOVE/samplesToMerge.list ${runDir}/fromSMOOVE/samples_merged_INV.raw.vcf.gz -Oz > ${runDir}/fromSMOOVE/samples_merged_INV.filt2.vcf.gz

# Filter variants with no ALT allele in any sample. Filter for call rate. Sort variants and clean header.
bcftools view -c1 ${runDir}/fromSMOOVE/samples_merged_DEL.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_DEL.filt3.vcf.gz
bcftools view -c1 ${runDir}/fromSMOOVE/samples_merged_DUP.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_DUP.filt3.vcf.gz 
bcftools view -c1 ${runDir}/fromSMOOVE/samples_merged_INV.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_INV.filt3.vcf.gz 
rm out.log

# Filter Variants by depth (Duphold)
${svpipeline_r_lib}/filterSVByFC.R ${runDir}/fromSMOOVE/samples_merged_DEL.filt3.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DEL.filt4.vcf.gz
${svpipeline_r_lib}/filterSVByFC.R ${runDir}/fromSMOOVE/samples_merged_DUP.filt3.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DUP.filt4.vcf.gz
${svpipeline_r_lib}/filterSVByFC.R ${runDir}/fromSMOOVE/samples_merged_INV.filt3.vcf.gz ${runDir}/fromSMOOVE/samples_merged_INV.filt4.vcf.gz

# Rename IDs
${svpipeline_r_lib}/rename_SVIds.R ${runDir}/fromSMOOVE/samples_merged_DEL.filt4.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DEL.filt5.vcf.gz
${svpipeline_r_lib}/rename_SVIds.R ${runDir}/fromSMOOVE/samples_merged_DUP.filt4.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DUP.filt5.vcf.gz
${svpipeline_r_lib}/rename_SVIds.R ${runDir}/fromSMOOVE/samples_merged_INV.filt4.vcf.gz ${runDir}/fromSMOOVE/samples_merged_INV.filt5.vcf.gz

# Remove some INFO tags 
zcat ${runDir}/fromSMOOVE/samples_merged_DEL.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > ${runDir}/fromSMOOVE/samples_merged_DEL.filt6.vcf
zcat ${runDir}/fromSMOOVE/samples_merged_DUP.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > ${runDir}/fromSMOOVE/samples_merged_DUP.filt6.vcf
zcat ${runDir}/fromSMOOVE/samples_merged_INV.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > ${runDir}/fromSMOOVE/samples_merged_INV.filt6.vcf

# Fix VCF header
java -jar $PICARD FixVcfHeader I=${runDir}/fromSMOOVE/samples_merged_DEL.filt6.vcf O=${runDir}/fromSMOOVE/samples_merged_DEL.filt7.vcf
java -jar $PICARD FixVcfHeader I=${runDir}/fromSMOOVE/samples_merged_DUP.filt6.vcf O=${runDir}/fromSMOOVE/samples_merged_DUP.filt7.vcf
java -jar $PICARD FixVcfHeader I=${runDir}/fromSMOOVE/samples_merged_INV.filt6.vcf O=${runDir}/fromSMOOVE/samples_merged_INV.filt7.vcf

# Compress results
parallel bgzip -f {} ::: ${runDir}/fromSMOOVE/*.vcf
rm ${runDir}/fromSMOOVE/*.idx

# Add variant info (Final Results)
bcftools +fill-tags -Oz -o ${runDir}/fromSMOOVE/samples_merged_DEL.Final.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DEL.filt7.vcf.gz
bcftools +fill-tags -Oz -o ${runDir}/fromSMOOVE/samples_merged_DUP.Final.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DUP.filt7.vcf.gz
bcftools +fill-tags -Oz -o ${runDir}/fromSMOOVE/samples_merged_INV.Final.vcf.gz ${runDir}/fromSMOOVE/samples_merged_INV.filt7.vcf.gz

# Filter by MAF (5%)
bcftools view -q 0.05:minor ${runDir}/fromSMOOVE/samples_merged_DEL.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_DEL.rate90.maf05.Final.vcf.gz
bcftools view -q 0.05:minor ${runDir}/fromSMOOVE/samples_merged_DUP.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_DUP.rate90.maf05.Final.vcf.gz
bcftools view -q 0.05:minor ${runDir}/fromSMOOVE/samples_merged_INV.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_INV.rate90.maf05.Final.vcf.gz

# Filter by MAF (1%)
bcftools view -q 0.01:minor ${runDir}/fromSMOOVE/samples_merged_DEL.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_DEL.rate90.maf01.Final.vcf.gz
bcftools view -q 0.01:minor ${runDir}/fromSMOOVE/samples_merged_DUP.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_DUP.rate90.maf01.Final.vcf.gz
bcftools view -q 0.01:minor ${runDir}/fromSMOOVE/samples_merged_INV.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o ${runDir}/fromSMOOVE/samples_merged_INV.rate90.maf01.Final.vcf.gz

# Compress vcf files
parallel bgzip -f -i {} ::: ${runDir}/fromSMOOVE/*.vcf
parallel tabix -f -p vcf {} ::: ${runDir}/fromSMOOVE/*.Final.vcf.gz

# Get tables for plots
bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' ${runDir}/fromSMOOVE/samples_merged_DEL.Final.vcf.gz > ${runDir}/fromSMOOVE/DEL.dat
bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' ${runDir}/fromSMOOVE/samples_merged_DUP.Final.vcf.gz > ${runDir}/fromSMOOVE/DUP.dat
bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' ${runDir}/fromSMOOVE/samples_merged_INV.Final.vcf.gz > ${runDir}/fromSMOOVE/INV.dat

sed -e 's/^/chr/' ${runDir}/fromSMOOVE/DEL.dat > ${runDir}/fromSMOOVE/DEL.chr.dat
bedtools annotate -both -i ${runDir}/fromSMOOVE/DEL.chr.dat -files ${sd} ${sr} > ${runDir}/fromSMOOVE/DEL.sd_sr_cov.txt
sed -e 's/^/chr/' ${runDir}/fromSMOOVE/DUP.dat > ${runDir}/fromSMOOVE/DUP.chr.dat
bedtools annotate -both -i ${runDir}/fromSMOOVE/DUP.chr.dat -files ${sd} ${sr} > ${runDir}/fromSMOOVE/DUP.sd_sr_cov.txt
sed -e 's/^/chr/' ${runDir}/fromSMOOVE/INV.dat > ${runDir}/fromSMOOVE/INV.chr.dat
bedtools annotate -both -i ${runDir}/fromSMOOVE/INV.chr.dat -files ${sd} ${sr} > ${runDir}/fromSMOOVE/INV.sd_sr_cov.txt

# Concat sv types (including MELT INS)
bcftools concat -a -Oz -o ${runDir}/fromSMOOVE/samples_merged_ALL.Final.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DEL.Final.vcf.gz ${runDir}/fromSMOOVE/samples_merged_DUP.Final.vcf.gz ${runDir}/fromSMOOVE/samples_merged_INV.Final.vcf.gz ${runDir}/fromMELT/samples_merged_INS.Final.vcf.gz ${runDir}/fromSURVIVOR/samples_merged_TRA_BND.Final.vcf.gz
tabix -p vcf ${runDir}/fromSMOOVE/samples_merged_ALL.Final.vcf.gz
