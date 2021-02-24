#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_SVE.sh
## Usage: ./run_02_SVE.sh myConfigFile.config

## This script will run BreakDancer, BreakSeq2 and CNVnator using SVE pipeline.
## NOTE: After finishing running all jobs in this step run collect_02_SVE.sh to collect the results.

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

## Run folder
mkdir -p ${outputFolderPath} || exit 1
mkdir -p ${outputFolderPath}/02_SVtools || exit 1
runDir=${outputFolderPath}"/02_SVtools/SVEPipeline"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/alignments || exit 1
mkdir -p ${runDir}/svCalls || exit 1
mkdir -p ${runDir}/logs || exit 1
cd ${runDir}

module purge
module load sve

# =============================================================================
#                    Loop over each sample and run SV discovery
# =============================================================================
## Call SVs (3 algorithms, BreakDancer, BreakSeq2 and CNVnator)
# Iterate over bam files list
while read bamFile; do

   name=$(basename "$bamFile")
   id=`echo $name | cut -d '.' -f 1`
   filename=$(basename -- "$bamFile")
   extension="${filename##*.}"
   filename="${filename%.*}"

   echo "|----- $bamFile ----- $id ----- $referenceFasta -----|"

   if [ ! -s "${runDir}/svCalls/${filename}_S4.vcf" ]
   then
   ## Breakdancer - OK
   echo "sve call -r ${referenceFasta} \
      -o ${runDir}/svCalls \
      -t 1 \
      -M 8 \
      -a breakdancer \
      -g hg19 ${bamFile}" \
   | bsub -g /cnv/sve \
      -n 1 \
      -R "rusage[mem=16384]" \
      -W 144:00 \
      -oo ${runDir}/logs/SVE.breakdancer.${id}.out \
      -eo ${runDir}/logs/SVE.breakdancer.${id}.err \
      -P acc_ad-omics \
      -q premium \
      -m bode \
      -J ${id}.SVE.breakdancer
   fi

   if [ ! -s "${runDir}/svCalls/${filename}_S10.vcf" ]
   then
   ## CNVnator - OK
   echo "sve call -r ${referenceFasta} \
      -o ${runDir}/svCalls \
      -t 1 \
      -M 8 \
      -a cnvnator \
      -g hg19 ${bamFile}" \
   | bsub -g /cnv/sve \
      -n 1 \
      -R "rusage[mem=16384]" \
      -W 144:00 \
      -oo ${runDir}/logs/SVE.cnvnator.${id}.out \
      -eo ${runDir}/logs/SVE.cnvnator.${id}.err \
      -P acc_ad-omics \
      -q premium \
      -m bode \
      -J ${id}.SVE.cnvnator
   fi

   if [ ! -s "${runDir}/svCalls/${filename}_S35.vcf" ]
   then
   ## Breakseq - OK with local python
   echo "sve call -r ${referenceFasta} \
      -o ${runDir}/svCalls \
      -t 1 \
      -M 8 \
      -a breakseq \
      -g hg19 ${bamFile}" \
   | bsub -g /cnv/sve \
      -n 1 \
      - "rusage[mem=16384]" \
      -W 144:00 \
      -oo ${runDir}/logs/SVE.breakseq.${id}.out \
      -eo ${runDir}/logs/SVE.breakseq.${id}.err \
      -P acc_ad-omics \
      -q premium \
      -m bode \
      -J ${id}.SVE.breakseq
   fi

done < $bamFileList


