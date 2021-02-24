#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 00_GeneticPC.sh
## Usage: ./collect_00_GeneticPC.sh myConfigFile.config

## This script must be ran after all jobs from run_00_GeneticPC.sh are finished.
## Final results will be stored inside the folder 00_GeneticPC.

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

# Run folder
mkdir -p ${outputFolderPath} || exit 1
runDir=${outputFolderPath}"/00_GeneticPC"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1
mkdir -p ${runDir}/extracted || exit 1

# Loading modules
ml purge
ml load htslib

cd ${runDir}

${path_to_somalier}/somalier relate ${runDir}/extracted/*.somalier

python ${path_to_somalier}/Repo/somalier/scripts/ancestry-predict.py --labels ${path_to_somalier}/Repo/somalier/scripts/ancestry-labels-1kg.tsv --samples ${runDir}/extracted/*.somalier --backgrounds ${path_to_somalier}/1kg-somalier/*.somalier --plot ${runDir}/pca.pdf

${path_to_somalier}/somalier ancestry pca --labels ${path_to_somalier}/Repo/somalier/scripts/ancestry-labels-1kg.tsv ${path_to_somalier}/1kg-somalier/*.somalier ${runDir}/extracted/*.somalier
