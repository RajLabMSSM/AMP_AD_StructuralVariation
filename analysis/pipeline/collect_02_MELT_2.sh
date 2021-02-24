#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_MELT
## Usage: ./collect_02_MELT_2.sh myConfigFile.config

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
#                           Genotype
# =============================================================================

while [ 1 ]
do

  read entry || break

  id=$(echo $entry | awk -F ' ' '{print $1}')
  bam=$(echo $entry | awk -F' ' '{print $2}')
  snv=$(echo $entry | awk -F' ' '{print $3}')
  sex=$(echo $entry | awk -F' ' '{print $4}')

  echo "|-- $bam ---- $id --|"

  sampleDir=${runDir}/${id}
  mkdir -p ${sampleDir}

  ##### Genotype ALU #####
  echo "java -Xmx2G -jar ${MELT_DIR}/MELT.jar Genotype \\" > ${sampleDir}/gt_alu_MELT.sh
  echo "  -bamfile ${bam} \\" >> ${sampleDir}/gt_alu_MELT.sh
  echo "  -t ${runDir}/ALU.txt \\" >> ${sampleDir}/gt_alu_MELT.sh
  echo "  -h ${referenceFasta} \\" >> ${sampleDir}/gt_alu_MELT.sh
  echo "  -w ${sampleDir} \\" >> ${sampleDir}/gt_alu_MELT.sh
  echo "  -p ${runDir}" >> ${sampleDir}/gt_alu_MELT.sh

  chmod a+x ${sampleDir}/gt_alu_MELT.sh

  if [ ! -e ${sampleDir}/${id}.final.ALU.tsv ]; then
    bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" -W 144:00 -oo ${sampleDir}/${id}.alu.out -eo ${sampleDir}/${id}.alu.err -P acc_ad-omics -q premium -J MELT_${id} ${sampleDir}/gt_alu_MELT.sh
  else
    if [ "$(wc -l < ${runDir}/ALU.pre_geno.tsv)" -eq "$(wc -l < ${sampleDir}/${id}.final.ALU.tsv)" ]; then 
      echo 'Sample already processed!'; 
    else 
      bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" -W 144:00 -oo ${sampleDir}/${id}.alu.out -eo ${sampleDir}/${id}.alu.err -P acc_ad-omics -q premium -J MELT_${id} ${sampleDir}/gt_alu_MELT.sh
    fi
  fi


  ##### Genotype SVA #####
  echo "java -Xmx2G -jar ${MELT_DIR}/MELT.jar Genotype \\" > ${sampleDir}/gt_sva_MELT.sh
  echo "  -bamfile ${bam} \\" >> ${sampleDir}/gt_sva_MELT.sh
  echo "  -t ${runDir}/SVA.txt \\" >> ${sampleDir}/gt_sva_MELT.sh
  echo "  -h ${referenceFasta} \\" >> ${sampleDir}/gt_sva_MELT.sh
  echo "  -w ${sampleDir} \\" >> ${sampleDir}/gt_sva_MELT.sh
  echo "  -p ${runDir}" >> ${sampleDir}/gt_sva_MELT.sh

  chmod a+x ${sampleDir}/gt_sva_MELT.sh

  if [ ! -e ${sampleDir}/${id}.final.SVA.tsv ]; then
    bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" -W 144:00 -oo ${sampleDir}/${id}.sva.out -eo ${sampleDir}/${id}.sva.err -P acc_ad-omics -q premium -J MELT_${id} ${sampleDir}/gt_sva_MELT.sh
  else
    if [ "$(wc -l < ${runDir}/SVA.pre_geno.tsv)" -eq "$(wc -l < ${sampleDir}/${id}.final.SVA.tsv)" ]; then 
      echo 'Sample already processed!'; 
    else 
      bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" -W 144:00 -oo ${sampleDir}/${id}.sva.out -eo ${sampleDir}/${id}.sva.err -P acc_ad-omics -q premium -J MELT_${id} ${sampleDir}/gt_sva_MELT.sh
    fi  
  fi

  ##### Genotype LINE1 #####
  echo "java -Xmx2G -jar ${MELT_DIR}/MELT.jar Genotype \\" > ${sampleDir}/gt_line1_MELT.sh
  echo "  -bamfile ${bam} \\" >> ${sampleDir}/gt_line1_MELT.sh
  echo "  -t ${runDir}/LINE1.txt \\" >> ${sampleDir}/gt_line1_MELT.sh
  echo "  -h ${referenceFasta} \\" >> ${sampleDir}/gt_line1_MELT.sh
  echo "  -w ${sampleDir} \\" >> ${sampleDir}/gt_line1_MELT.sh
  echo "  -p ${runDir}" >> ${sampleDir}/gt_line1_MELT.sh

  chmod a+x ${sampleDir}/gt_line1_MELT.sh

  if [ ! -e ${sampleDir}/${id}.final.LINE1.tsv ]; then
    bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" -W 144:00 -oo ${sampleDir}/${id}.line1.out -eo ${sampleDir}/${id}.line1.err -P acc_ad-omics -q premium -J MELT_${id} ${sampleDir}/gt_line1_MELT.sh
  else
    if [ "$(wc -l < ${runDir}/LINE1.pre_geno.tsv)" -eq "$(wc -l < ${sampleDir}/${id}.final.LINE1.tsv)" ]; then 
      echo 'Sample already processed!'; 
    else
      bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=2000]" -W 144:00 -oo ${sampleDir}/${id}.line1.out -eo ${sampleDir}/${id}.line1.err -P acc_ad-omics -q premium -J MELT_${id} ${sampleDir}/gt_line1_MELT.sh
    fi      
  fi

done < $id_mapping

