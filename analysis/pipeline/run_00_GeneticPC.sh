#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 00_GeneticPC.sh
## Usage: ./run_00_GeneticPC.sh myConfigFile.config

## Script to get sample relatedness and genetic ancestry directly from BAM files using somalier.
## A folder named "00_GeneticPC" will be created

## NOTE: After finishing running all jobs in this step run collect_00_GeneticPC.sh to collect the results.

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

# Run folder
mkdir -p ${outputFolderPath} || exit 1
runDir=${outputFolderPath}"/00_GeneticPC"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1
mkdir -p ${runDir}/extracted || exit 1

# Loading modules
ml purge
ml load htslib

# =============================================================================
#                    Loop over each sample to extract sites
# =============================================================================

while [ 1 ]
do
    read line || break
    name=$(basename "$line")
    id=`echo $name | cut -d '.' -f 1`

    echo "${path_to_somalier}/somalier extract -d ${runDir}/extracted/ --sites ${path_to_somalier}/sites.GRCh37.vcf.gz -f ${referenceFasta} ${line}" | bsub -n 1 -R "rusage[mem=4000]" -W 12:00 -oo ${runDir}/logs/${id}.out -eo ${runDir}/logs/${id}.err -P acc_ad-omics -q express -J extract.${id}

done < $bamFileList
