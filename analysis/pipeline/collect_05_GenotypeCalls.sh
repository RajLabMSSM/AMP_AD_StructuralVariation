#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 05_GenotypeCalls
## Usage: ./collect_05_GenotypeCalls.sh myConfigFile.config

## This script must be ran after all jobs from run_05_GenotypeCalls.sh are finished.

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
# Collect Smoove results and prepare vcf files with final genotypes
##

prevResultsDir=${outputFolderPath}"/04_GenotypedCalls_Smoove"
runDir=${prevResultsDir}"/Results"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/gnomadSVs || exit 1
mkdir -p ${runDir}/mergedCohorts || exit 1
mkdir -p ${runDir}/denovoSVs || exit 1

ml purge
ml R/3.6.2

# =============================================================================
#                      Merged genotyped results
# =============================================================================

# Check if all samples were completed
del_lines=$(grep -v "^#" ${prevResultsDir}/samples_merged_DEL.vcf | wc -l | cut -f1 -d' ')
dup_lines=$(grep -v "^#" ${prevResultsDir}/samples_merged_DUP.vcf | wc -l | cut -f1 -d' ')
inv_lines=$(grep -v "^#" ${prevResultsDir}/samples_merged_INV.vcf | wc -l | cut -f1 -d' ')
while [ 1 ]
do

    read entry || break

    id=$(echo $entry | awk -F ' ' '{print $1}')
    bam=$(echo $entry | awk -F' ' '{print $2}')
    snv=$(echo $entry | awk -F' ' '{print $3}')
    sex=$(echo $entry | awk -F' ' '{print $4}')

    sample_del_lines=$(zcat ${prevResultsDir}/DELs/${id}/${id}-smoove.genotyped.vcf.gz | grep -v "^#" | wc -l | cut -f1 -d' ')
    if [ ! "$del_lines" -eq "$sample_del_lines" ]; then
        echo "WARNING SAMPLE FAILED -- DELs: $id"
    fi
    sample_dup_lines=$(zcat ${prevResultsDir}/DUPs/${id}/${id}-smoove.genotyped.vcf.gz | grep -v "^#" | wc -l | cut -f1 -d' ')
    if [ ! "$dup_lines" -eq "$sample_dup_lines" ]; then
        echo "WARNING SAMPLE FAILED -- DUPs: $id"
    fi
    sample_inv_lines=$(zcat ${prevResultsDir}/INVs/${id}/${id}-smoove.genotyped.vcf.gz | grep -v "^#" | wc -l | cut -f1 -d' ')
    if [ ! "$inv_lines" -eq "$sample_inv_lines" ]; then
        echo "WARNING SAMPLE FAILED -- INVs: $id"
    fi

done < $id_mapping

# paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.
echo "smoove paste --name DEL --outdir ${runDir} ${prevResultsDir}/DELs/*/*.genotyped.vcf.gz" | bsub -n 1 -R "rusage[mem=64000]" -W 12:00 -oo ${runDir}/paste.del.out -eo ${runDir}/paste.del.err -P acc_ad-omics -q express -J merge_del
echo "smoove paste --name DUP --outdir ${runDir} ${prevResultsDir}/DUPs/*/*.genotyped.vcf.gz" | bsub -n 1 -R "rusage[mem=64000]" -W 12:00 -oo ${runDir}/paste.dup.out -eo ${runDir}/paste.dup.err -P acc_ad-omics -q express -J merge_dup
echo "smoove paste --name INV --outdir ${runDir} ${prevResultsDir}/INVs/*/*.genotyped.vcf.gz" | bsub -n 1 -R "rusage[mem=64000]" -W 12:00 -oo ${runDir}/paste.inv.out -eo ${runDir}/paste.inv.err -P acc_ad-omics -q express -J merge_inv
