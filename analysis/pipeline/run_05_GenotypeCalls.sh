#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 05_GenotypeCalls
## Usage: ./run_05_GenotypeCalls.sh myConfigFile.config

## Script to get genotype population merged calls using smoove. 
## A job is submitted per sample. Only DEL, DUP and INV are genotyped using this method.
## A folder named "05_GenotypeCalls" will be created

## NOTE: After finishing running all jobs in this step run collect_05_GenotypeCalls.sh to collect the results.

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
# This script will call genotypes using the tool smoove (SVtyper).

ml purge
ml R/3.6.2 

## Run folder
mkdir -p ${outputFolderPath} || exit 1
runDir=${outputFolderPath}"/05_GenotypeCalls"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1

# =============================================================================
#                                   Prepare data
# =============================================================================

## Folder with SVs for each sample (after merging tools with SURVIVOR)
callsDir=${outputFolderPath}/04_MergeSamples/population

## Filter per size (>50 bp AND <50M bp)
SURVIVOR filter ${callsDir}/samples_merged_ALL.vcf NA 50 50000000 0 -1 ${runDir}/samples_merged_ALL.filt.vcf

## Split VCF into DEL, DUP, INV and cut out individual information
grep -E '^#|SVTYPE=DEL' ${runDir}/samples_merged_ALL.filt.vcf | cut -f1,2,3,4,5,6,7,8,9 > ${runDir}/samples_merged_DEL.vcf
grep -E '^#|SVTYPE=DUP' ${runDir}/samples_merged_ALL.filt.vcf | cut -f1,2,3,4,5,6,7,8,9 > ${runDir}/samples_merged_DUP.vcf
grep -E '^#|SVTYPE=INV' ${runDir}/samples_merged_ALL.filt.vcf | cut -f1,2,3,4,5,6,7,8,9 > ${runDir}/samples_merged_INV.vcf

## Add required CIPOS95 and CIEND95 fields
# DEL
${svpipeline_r_lib}/fixCIPOS2.R ${runDir}/samples_merged_DEL.vcf ${runDir}/samples_merged_DEL.vcf.gz
gunzip -f ${runDir}/samples_merged_DEL.vcf.gz
mkdir -p ${runDir}/DELs

# DUP
${svpipeline_r_lib}/fixCIPOS2.R ${runDir}/samples_merged_DUP.vcf ${runDir}/samples_merged_DUP.vcf.gz
gunzip -f ${runDir}/samples_merged_DUP.vcf.gz
mkdir -p ${runDir}/DUPs

# INV
sed -i 's/<BND>/<INV>/g' ${runDir}/samples_merged_INV.vcf
${svpipeline_r_lib}/fixCIPOS2.R ${runDir}/samples_merged_INV.vcf ${runDir}/samples_merged_INV.vcf.gz
gunzip -f ${runDir}/samples_merged_INV.vcf.gz
mkdir -p ${runDir}/INVs

# Get number of SVs to be genotyped
del_lines=$(grep -v "^#" ${runDir}/samples_merged_DEL.vcf | wc -l | cut -f1 -d' ')
dup_lines=$(grep -v "^#" ${runDir}/samples_merged_DUP.vcf | wc -l | cut -f1 -d' ')
inv_lines=$(grep -v "^#" ${runDir}/samples_merged_INV.vcf | wc -l | cut -f1 -d' ')
cat /dev/null > ${runDir}/samples.report

# =============================================================================
#                    Loop over each sample to genotype
# =============================================================================

while [ 1 ]
do

    read entry || break

    id=$(echo $entry | awk -F ' ' '{print $1}')
    bam=$(echo $entry | awk -F' ' '{print $2}')
    snv=$(echo $entry | awk -F' ' '{print $3}')
    sex=$(echo $entry | awk -F' ' '{print $4}')

    echo $id

    sample_del_lines=$(zcat ${runDir}/DELs/${id}/${id}-smoove.genotyped.vcf.gz | grep -v "^#" | wc -l | cut -f1 -d' ')
    if [ ! "$del_lines" -eq "$sample_del_lines" ]; then
        echo "DELs: $id" >> ${runDir}/samples.report
        echo "smoove genotype -p 1 --name $id --outdir ${runDir}/DELs/$id --fasta $referenceFasta --duphold --vcf ${runDir}/samples_merged_DEL.vcf $bam" | bsub -g /smoove -n 1 -R "rusage[mem=4000]" -W 12:00 -oo ${runDir}/logs/${id}.del.out -eo ${runDir}/logs/${id}.del.err -P acc_ad-omics -q express -J del_${id}
    fi
    sample_dup_lines=$(zcat ${runDir}/DUPs/${id}/${id}-smoove.genotyped.vcf.gz | grep -v "^#" | wc -l | cut -f1 -d' ')
    if [ ! "$dup_lines" -eq "$sample_dup_lines" ]; then
        echo "DUPs: $id" >> ${runDir}/samples.report
        echo "smoove genotype -p 1 --name $id --outdir ${runDir}/DUPs/$id --fasta $referenceFasta --duphold --vcf ${runDir}/samples_merged_DUP.vcf $bam" | bsub -g /smoove -n 1 -R "rusage[mem=4000]" -W 12:00 -oo ${runDir}/logs/${id}.dup.out -eo ${runDir}/logs/${id}.dup.err -P acc_ad-omics -q express -J dup_${id}
    fi
    sample_inv_lines=$(zcat ${runDir}/INVs/${id}/${id}-smoove.genotyped.vcf.gz | grep -v "^#" | wc -l | cut -f1 -d' ')
    if [ ! "$inv_lines" -eq "$sample_inv_lines" ]; then
        echo "INVs: $id" >> ${runDir}/samples.report
        echo "smoove genotype -p 1 --name $id --outdir ${runDir}/INVs/$id --fasta $referenceFasta --duphold --vcf ${runDir}/samples_merged_INV.vcf $bam" | bsub -g /smoove -n 1 -R "rusage[mem=4000]" -W 12:00 -oo ${runDir}/logs/${id}.inv.out -eo ${runDir}/logs/${id}.inv.err -P acc_ad-omics -q express -J inv_${id}
    fi

done < $id_mapping
