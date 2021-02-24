#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_Delly.sh
## Usage: ./run_02_Delly.sh myConfigFile.config

## NOTE: After finishing running all jobs in this step run collect_02_Delly.sh to collect the results.

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
runDir=${outputFolderPath}"/02_SVtools/Delly"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/data || exit 1
mkdir -p ${runDir}/logs || exit 1

module purge
module load bcftools
module load vcftools
module load anaconda3
module load delly

# =============================================================================
#                    Loop over each sample and run SV discovery
# =============================================================================

while [ 1 ]
do

   read bamFile || break
   fileName=$(basename "$bamFile")
   id=`echo $fileName | cut -d '.' -f 1`
   echo "|-- $bamFile ---- $id --|"

   if [ ! -s "${runDir}/data/${id}.DEL.bcf" ]; then
      echo "delly call -g ${referenceFasta} \
         -o ${runDir}/data/${id}.DEL.bcf \
         -x ${delly_exclude} \
         -t DEL -e -i ${bamFile}" \
      | bsub -g /cnv/delly \
         -q premium \
         -P acc_ad-omics \
         -oo ${runDir}/logs/${id}.DEL.out \
         -eo ${runDir}/logs/${id}.DEL.err \
         -W 72:00 \
         -n 1 \
         -R "rusage[mem=12000]"
   fi
   if [ ! -s "${runDir}/data/${id}.INS.bcf" ]; then
      echo "delly call -g ${referenceFasta} \
         -o ${runDir}/data/${id}.INS.bcf \
         -x ${delly_exclude} \
         -t INS -e -i ${bamFile}" \
      | bsub -g /cnv/delly \
         -q premium \
         -P acc_ad-omics \
         -oo ${runDir}/logs/${id}.INS.out \
         -eo ${runDir}/logs/${id}.INS.err \
         -W 72:00 \
         -n 1 \
         -R "rusage[mem=12000]"
   fi
   if [ ! -s "${runDir}/data/${id}.DUP.bcf" ]; then
      echo "delly call -g ${referenceFasta} \
         -o ${runDir}/data/${id}.DUP.bcf \
         -x ${delly_exclude} \
         -t DUP -e -i ${bamFile}" \
      | bsub -g /cnv/delly \
         -q premium \
         -P acc_ad-omics \
         -oo ${runDir}/logs/${id}.DUP.out \
         -eo ${runDir}/logs/${id}.DUP.err \
         -W 72:00 \
         -n 1 \
         -R "rusage[mem=12000]"
   fi
   if [ ! -s "${runDir}/data/${id}.INV.bcf" ]; then
      echo "delly call -g ${referenceFasta} \
         -o ${runDir}/data/${id}.INV.bcf \
         -x ${delly_exclude} \
         -t INV -e -i ${bamFile}" \
      | bsub -g /cnv/delly \
         -q premium \
         -P acc_ad-omics \
         -oo ${runDir}/logs/${id}.INV.out \
         -eo ${runDir}/logs/${id}.INV.err \
         -W 72:00 \
         -n 1 \
         -R "rusage[mem=12000]"
   fi
   if [ ! -s "${runDir}/data/${id}.BND.bcf" ]; then
      echo "delly call -g ${referenceFasta} \
         -o ${runDir}/data/${id}.BND.bcf \
         -x ${delly_exclude} \
         -t BND -e -i ${bamFile}" \
      | bsub -g /cnv/delly \
         -q premium \
         -P acc_ad-omics \
         -oo ${runDir}/logs/${id}.BND.out \
         -eo ${runDir}/logs/${id}.BND.err \
         -W 72:00 \
         -n 1 \
         -R "rusage[mem=12000]"
   fi

done < $bamFileList
