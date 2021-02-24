#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_MELT
## Usage: ./collect_02_MELT_3.sh myConfigFile.config

## This script must be ran after all jobs from run_02_MELT.sh are finished.
## A folder Results will be created inside the folder 02_SVtools/MELT.

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

mkdir -p ${outputFolderPath} || exit 1
mkdir -p ${outputFolderPath}/02_SVtools || exit 1
prevDir=${outputFolderPath}"/02_SVtools/MELT" 
mkdir -p ${prevDir} || exit 1
runDir=${prevDir}"/Results"
mkdir -p ${runDir} || exit 1

## Required modules. 
module purge
module load java
module load bowtie2

MELT_DIR="/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5"

echo "/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5/me_refs/1KGP_Hg19/ALU_MELT.zip" > ${runDir}/ALU.txt
echo "/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5/me_refs/1KGP_Hg19/SVA_MELT.zip" > ${runDir}/SVA.txt
echo "/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip" > ${runDir}/LINE1.txt

# =============================================================================
#                           MakeVCF
# =============================================================================

## ALU
java -Xmx2G -jar ${MELT_DIR}/MELT.jar MakeVCF \
    -genotypingdir ${runDir}/ \
    -h ${referenceFasta} \
    -t ${runDir}/ALU.txt \
    -w ${runDir}/ \
    -p ${runDir}/ \
    -o ${runDir}

## SVA
java -Xmx2G -jar ${MELT_DIR}/MELT.jar MakeVCF \
    -genotypingdir ${runDir}/ \
    -h ${referenceFasta} \
    -t ${runDir}/SVA.txt \
    -w ${runDir}/ \
    -p ${runDir}/ \
    -o ${runDir}

## LINE1
java -Xmx2G -jar ${MELT_DIR}/MELT.jar MakeVCF \
    -genotypingdir ${runDir}/ \
    -h ${referenceFasta} \
    -t ${runDir}/LINE1.txt \
    -w ${runDir}/ \
    -p ${runDir}/ \
    -o ${runDir}
