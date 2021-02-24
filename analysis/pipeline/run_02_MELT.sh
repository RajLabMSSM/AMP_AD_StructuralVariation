#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_MELT.sh
## Usage: ./run_02_MELT.sh myConfigFile.config

## NOTE: After finishing running all jobs in this step run collect_02_MELT_1.sh, collect_02_MELT_2.sh, and collect_02_MELT_3.sh to collect the results.

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
runDir=${outputFolderPath}"/02_SVtools/MELT"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1

## Required modules.
module purge
module load java
module load bowtie2

MELT_DIR="/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5"

echo "/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5/me_refs/1KGP_Hg19/ALU_MELT.zip" > ${runDir}/mei_list.txt
echo "/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5/me_refs/1KGP_Hg19/SVA_MELT.zip" >> ${runDir}/mei_list.txt
echo "/hpc/users/viallr01/ad-omics/ricardo/MyApps/MELT/MELTv2.1.5/me_refs/1KGP_Hg19/LINE1_MELT.zip" >> ${runDir}/mei_list.txt

# =============================================================================
#                    Loop over each sample and run SV discovery
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
  mkdir -p ${sampleDir}/ALU
  mkdir -p ${sampleDir}/LINE1
  mkdir -p ${sampleDir}/SVA

  cat /dev/null > ${sampleDir}/runMELT.sh

  echo "## Preprocessing .bam Files for MELT" >> ${sampleDir}/runMELT.sh
  echo "java -Xmx2G -jar ${MELT_DIR}/MELT.jar Preprocess -bamfile ${bam} -h $referenceFasta" >> ${sampleDir}/runMELT.sh
  echo "" >> ${sampleDir}/runMELT.sh
  echo "## Running MELT Using MELT-SPLIT" >> ${sampleDir}/runMELT.sh
  echo "## IndivAnalysis – MEI discovery in individual samples" >> ${sampleDir}/runMELT.sh
  echo "## ALU" >> ${sampleDir}/runMELT.sh
  echo "java -Xmx6G -jar ${MELT_DIR}/MELT.jar IndivAnalysis \\" >> ${sampleDir}/runMELT.sh
  echo "  -c 30 \\" >> ${sampleDir}/runMELT.sh
  echo "  -h ${referenceFasta} \\" >> ${sampleDir}/runMELT.sh
  echo "  -bamfile ${bam} \\" >> ${sampleDir}/runMELT.sh
  echo "  -t ${MELT_DIR}/me_refs/1KGP_Hg19/ALU_MELT.zip \\" >> ${sampleDir}/runMELT.sh
  echo "  -w ${sampleDir}/ALU/" >> ${sampleDir}/runMELT.sh
  echo "rm -rf ${sampleDir}/ALU/*tmp" >> ${sampleDir}/runMELT.sh
  #echo "rm -rf ${sampleDir}/ALU/*bam" >> ${sampleDir}/runMELT.sh
  echo "" >> ${sampleDir}/runMELT.sh
  echo "## LINE1" >> ${sampleDir}/runMELT.sh
  echo "java -Xmx6G -jar ${MELT_DIR}/MELT.jar IndivAnalysis \\" >> ${sampleDir}/runMELT.sh
  echo "  -c 30 \\" >> ${sampleDir}/runMELT.sh
  echo "  -h ${referenceFasta} \\" >> ${sampleDir}/runMELT.sh
  echo "  -bamfile ${bam} \\" >> ${sampleDir}/runMELT.sh
  echo "  -t ${MELT_DIR}/me_refs/1KGP_Hg19/LINE1_MELT.zip \\" >> ${sampleDir}/runMELT.sh
  echo "  -w ${sampleDir}/LINE1/" >> ${sampleDir}/runMELT.sh
  echo "rm -rf ${sampleDir}/LINE1/*tmp" >> ${sampleDir}/runMELT.sh
  #echo "rm -rf ${sampleDir}/LINE1/*bam" >> ${sampleDir}/runMELT.sh
  echo "" >> ${sampleDir}/runMELT.sh
  echo "## SVA" >> ${sampleDir}/runMELT.sh
  echo "java -Xmx6G -jar ${MELT_DIR}/MELT.jar IndivAnalysis \\" >> ${sampleDir}/runMELT.sh
  echo "  -c 30 \\" >> ${sampleDir}/runMELT.sh
  echo "  -h ${referenceFasta} \\" >> ${sampleDir}/runMELT.sh
  echo "  -bamfile ${bam} \\" >> ${sampleDir}/runMELT.sh
  echo "  -t ${MELT_DIR}/me_refs/1KGP_Hg19/SVA_MELT.zip \\" >> ${sampleDir}/runMELT.sh
  echo "  -w ${sampleDir}/SVA/" >> ${sampleDir}/runMELT.sh
  echo "rm -rf ${sampleDir}/SVA/*tmp" >> ${sampleDir}/runMELT.sh
  #echo "rm -rf ${sampleDir}/SVA/*bam" >> ${sampleDir}/runMELT.sh
  echo "" >> ${sampleDir}/runMELT.sh
  echo "## Remove Preprocess files" >> ${sampleDir}/runMELT.sh
  echo "rm -rf ${bam}.d*" >> ${sampleDir}/runMELT.sh
  echo "rm -rf ${bam}.f*" >> ${sampleDir}/runMELT.sh

  chmod a+x ${sampleDir}/runMELT.sh

  if [ ! -e ${sampleDir}/SVA/${id}.SVA.aligned.final.sorted.bam ]; then
    bsub -g /melt -n 1 -R "span[hosts=1]" -R "rusage[mem=32000]" -W 12:00 -oo ${sampleDir}/${id}.out -eo ${sampleDir}/${id}.err -P acc_ad-omics -q express -J MELT_${id} ${sampleDir}/runMELT.sh
  else
    echo "OK"
  fi

done < $id_mapping
