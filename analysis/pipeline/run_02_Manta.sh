#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_Manta.sh
## Usage: ./run_02_Manta.sh myConfigFile.config

## NOTE: After finishing running all jobs in this step run collect_02_Manta.sh to collect the results.

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
mkdir -p ${outputFolderPath}/02_SVtools || exit 1
runDir=${outputFolderPath}"/02_SVtools/Manta"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/data || exit 1

module purge
module load manta

# =============================================================================
#                    Loop over each sample and run SV discovery
# =============================================================================

while [ 1 ]
do

   read bamFile || break
   fileName=$(basename "$bamFile")
   id=`echo $fileName | cut -d '.' -f 1`
   echo "|-- $bamFile ---- $id --|"

   if [ ! -s "${runDir}/data/${id}/results/variants/diploidSV.vcf.gz" ]; then
      rm -rf ${runDir}/data/${id}/
      configManta.py --bam ${bamFile} \
         --referenceFasta ${referenceFasta} \
         --callRegions ${intervalListBed} \
         --callMemMb 2000 \
         --runDir ${runDir}/data/${id}/
      echo "${runDir}/data/${id}/runWorkflow.py -m local -j 1 -g 4" > ${runDir}/data/${id}/submit_job.sh
      echo "rm -rf ${runDir}/data/${id}/workspace" >> ${runDir}/data/${id}/submit_job.sh
      chmod a+x ${runDir}/data/${id}/submit_job.sh

      bsub -g /cnv/manta -q express -P acc_ad-omics -oo ${runDir}/data/${id}/${id}.out -eo ${runDir}/data/${id}/${id}.err -W 12:00 -n 1 -R "rusage[mem=6000]" ${runDir}/data/${id}/submit_job.sh

   fi

done < $bamFileList
