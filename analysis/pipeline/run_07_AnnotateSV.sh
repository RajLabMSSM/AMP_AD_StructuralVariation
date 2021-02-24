#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 07_AnnotateSV
## Usage: ./run_07_AnnotateSV.sh myConfigFile.config

## A folder named "07_AnnotateSV" will be created

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
runDir=${outputFolderPath}/07_AnnotateSV
mkdir -p ${runDir} || exit 1

ml purge
ml samtools
ml bcftools
ml R/3.6.2

export ANNOTSV="/hpc/users/viallr01/ad-omics/ricardo/MyApps/AnnotSV"
export PATH="/hpc/users/viallr01/ad-omics/ricardo/MyApps/AnnotSV/bin/AnnotSV:$PATH"

# =============================================================================
#                           Annotate using AnnotSV
# =============================================================================
### SVs from Smoove

bgzip -d -c ${outputFolderPath}/06_Filtering/fromSMOOVE/samples_merged_ALL.Final.vcf.gz > ${runDir}/samples_merged_ALL.Final.vcf

SURVIVOR vcftobed ${runDir}/samples_merged_ALL.Final.vcf -99999999 99999999 ${runDir}/samples_merged_ALL.Final.bed
cat ${runDir}/samples_merged_ALL.Final.bed | cut -f1,2,5,7,11 | awk -F '\t' '{ $3 = ($3 == "0" ? $2+1 : $3) } 1' OFS='\t' > ${runDir}/samples_merged_ALL.Final.clean.bed

AnnotSV.tcl -SVinputFile ${runDir}/samples_merged_ALL.Final.clean.bed -outputDir ${runDir}/AnnotSV -SVinputInfo 1 -reciprocal yes -svtBEDcol 5

rm ${runDir}/samples_merged_ALL.Final.vcf ${runDir}/samples_merged_ALL.Final.bed ${runDir}/samples_merged_ALL.Final.clean.bed

# =============================================================================
#                         Complement annotation
# =============================================================================
${svpipeline_r_lib}/complementAnnotation.R ${runDir}/AnnotSV/samples_merged_ALL.Final.clean.annotated.tsv ${runDir}/AnnotSV/samples_merged_ALL.Final.clean.annotated.plus.tsv h37
