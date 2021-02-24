#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 01_SampleQC.sh
## Usage: ./run_01_SampleQC.sh myConfigFile.config

## Script to get QC metrics from BAM files
## For each sample, 7 jobs are submitted
## Runs the following:
## 1 - Picard CollectMultipleMetrics - output: ${sampleID}.metrics.insert_size_metrics
## 2 - Samtools flagstat - output: ${sampleID}.flagstat
## 3 - Bamtools stats - output: ${sampleID}.stats
## 4 - Picard CollectWgsMetrics - output: ${sampleID}.wgs
## 5 - Picard EstimateLibraryComplexity - output: ${sampleID}.complexity
## 6 - Sex Check - output: ${sampleID}.sexCheck
## 7 - Aneuploidy Check - output: ${sampleID}.aneuploidyCheck

## This step is based on the the first module from the HOLMES pipeline. Some of its scripts are included in the QC_Scripts folder.
## Collins RL, et al. Defining the spectrum of large inversions, complex structural variation, and chromothripsis in the morbid genome. Genome Biol. (2017)

## NOTE: After finishing running all jobs in this step run collect_01_SampleQC.sh to collect the results.

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
runDir=${outputFolderPath}"/01_SampleQC"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1

# Loading modules 
ml purge
ml sambamba
ml picard
ml samtools
ml bamtools

# =============================================================================
#                    Loop over each sample to collect metrics
# =============================================================================

while [ 1 ]
do
    read line || break
    name=$(basename "$line")
    id=`echo $name | cut -d '.' -f 1`

    # Submit library metrics QC jobs
    echo "#----------------------------------------------------------#"
    echo "Submit library metrics QC jobs"
    echo "#----------------------------------------------------------#"
    echo "$line $id $referenceFasta"

    if [ ! -s "${runDir}/${id}.metrics.insert_size_metrics" ]; then
        echo "java -Xmx4g -jar ${PICARD} CollectMultipleMetrics TMP_DIR=${tmpDir} I=${line} AS=true R=${referenceFasta} VALIDATION_STRINGENCY=SILENT PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics OUTPUT=${runDir}/${id}.metrics" | bsub -n 1 -R "rusage[mem=12000]" -W 48:00 -oo ${runDir}/logs/QC1.${id}.out -eo ${runDir}/logs/QC1.${id}.err -P acc_ad-omics -q premium -J QC1 -g /cnv/qc
    fi

    if [ ! -s "${runDir}/${id}.flagstat" ]; then
        echo "sambamba view -h -f bam -F 'not secondary_alignment' ${line} | samtools flagstat /dev/stdin > ${runDir}/${id}.flagstat" | bsub -n 1 -R "rusage[mem=12000]" -W 48:00 -oo ${runDir}/logs/QC2.${id}.out -eo ${runDir}/logs/QC2.${id}.err -P acc_ad-omics -q premium -J QC2 -g /cnv/qc
    fi

    if [ ! -s "${runDir}/${id}.stats" ]; then
        echo "sambamba view -h -f bam -F 'not secondary_alignment' ${line} | bamtools stats -in /dev/stdin -insert > ${runDir}/${id}.stats" | bsub -n 1 -R "rusage[mem=12000]" -W 48:00 -oo ${runDir}/logs/QC3.${id}.out -eo ${runDir}/logs/QC3.${id}.err -P acc_ad-omics -q premium -J QC3 -g /cnv/qc
    fi

    if [ ! -s "${runDir}/${id}.wgs" ]; then
        echo "java -Xmx4g -jar ${PICARD} CollectWgsMetrics TMP_DIR=${tmpDir} I=${line} O=${runDir}/${id}.wgs R=${referenceFasta} VALIDATION_STRINGENCY=SILENT" | bsub -n 1 -R "rusage[mem=12000]" -W 48:00 -oo ${runDir}/logs/QC4.${id}.out -eo ${runDir}/logs/QC4/${id}.err -P acc_ad-omics -q premium -J QC4 -g /cnv/qc
    fi

    if [ ! -s "${runDir}/${id}.complexity" ]; then
        echo "java -Xmx500g -jar ${PICARD} EstimateLibraryComplexity TMP_DIR=${tmpDir} I=${line} O=${runDir}/${id}.complexity VALIDATION_STRINGENCY=SILENT" | bsub -n 1 -R "rusage[mem=450000]" -R himem -W 48:00 -oo ${runDir}/logs/QC5.${id}.out -eo ${runDir}/logs/QC5.${id}.err -P acc_ad-omics -q premium -J QC5 -g /cnv/qc
    fi

    # Submit sex check
    echo "#----------------------------------------------------------#"
    echo "Submit sex check"
    echo "#----------------------------------------------------------#"

    if [ ! -s "${runDir}/${id}.sexCheck" ]; then
        # $1=ID, $2=bam, $3=DICT, $4=runDir
        echo "${sexCheck} ${id} ${line} ${DICT} ${runDir}" | bsub -n 1 -R "rusage[mem=12000]" -W 48:00 -oo ${runDir}/logs/QC6.${id}.out -eo ${runDir}/logs/QC6.${id}.err -P acc_ad-omics -q premium -J QC6 -g /cnv/qc
    fi

    # Submit aneuploidy check
    echo "#----------------------------------------------------------#"
    echo "Submit aneuploidy check"
    echo "#----------------------------------------------------------#"

    if [ ! -s "${runDir}/${id}.aneuploidyCheck" ]; then
        # $1=ID, $2=bam, $3=DICT, $4=runDir, $5=NMASK
        echo "${chrCopyCount} ${id} ${line} ${DICT} ${runDir} ${NMASK}" | bsub -n 1 -R "rusage[mem=12000]" -W 48:00 -oo ${runDir}/logs/QC7.${id}.out -eo ${runDir}/logs/QC7.${id}.err -P acc_ad-omics -q premium -J QC7 -g /cnv/qc
    fi

done < $bamFileList
