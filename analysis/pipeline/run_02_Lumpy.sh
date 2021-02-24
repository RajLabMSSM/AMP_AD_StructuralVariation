#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_Lumpy.sh
## Usage: ./run_02_Lumpy.sh myConfigFile.config

## NOTE: After finishing running all jobs in this step run collect_02_Lumpy.sh to collect the results.

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
runDir=${outputFolderPath}"/02_SVtools/Lumpy"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/logs || exit 1
mkdir -p ${runDir}/data || exit 1
mkdir -p ${runDir}/temp || exit 1

## Load modules
module purge # Clear current modules
module load lumpy # includes lumpy-express
module load speedseq
module load zlib

# =============================================================================
#                    Loop over each sample and run SV discovery
# =============================================================================

cd ${runDir}

## Prepare BAM files and run Lumpyexpress (Automated breakpoint detection for standard analyses) for each sample
while [ 1 ]
do

    read line || break
    name=$(basename "$line")    
    id=`echo $name | cut -d '.' -f 1`

    echo "|--- $line --- $id --- $referenceFasta ---|"

    if [[ ! -s "${runDir}/data/${id}.gt.vcf" ]]; then # If file is empty or does not exist
    
       ## Extract the discordant paired-end alignments.
       echo "samtools view -b -F 1294 $line > ${runDir}/data/${id}.discordants.unsorted.bam" > ${runDir}/logs/${id}-lumpy.sh

       ## Extract the split-read alignments
       echo "samtools view -h $line | {path_to_lumpy}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${runDir}/data/${id}.splitters.unsorted.bam" >> ${runDir}/logs/${id}-lumpy.sh

       ## Sort both alignments
       echo "samtools sort ${runDir}/data/${id}.discordants.unsorted.bam -o ${runDir}/data/${id}.discordants.bam -T ${id}" >> ${runDir}/logs/${id}-lumpy.sh
       echo "samtools sort ${runDir}/data/${id}.splitters.unsorted.bam -o ${runDir}/data/${id}.splitters.bam -T ${id}" >> ${runDir}/logs/${id}-lumpy.sh

       rnd=$(hexdump -n 16 -v -e '/1 "%02X"' /dev/urandom)

       echo "${lumpyexpress_path}/lumpyexpress -P -k -v -x ${speedseq_annotations}/ceph18.b37.lumpy.exclude.2014-01-15.bed -B $line -S ${runDir}/data/${id}.splitters.bam -D ${runDir}/data/${id}.discordants.bam -o ${runDir}/data/${id}.vcf -K ${speedseq_annotations}/lumpyexpress.config -T ${runDir}/temp/${id}_${rnd}/" >> ${runDir}/logs/${id}-lumpy.sh
    
       echo "svtyper-sso -i ${runDir}/data/${id}.vcf -o ${runDir}/data/${id}.gt.vcf -B $line -T $referenceFasta -l ${runDir}/data/${id}.gt.json --cores 1" >> ${runDir}/logs/${id}-lumpy.sh
    
       chmod a+x ${runDir}/logs/${id}-lumpy.sh

       bsub -g /cnv/lumpy -q express -P acc_ad-omics -oo ${runDir}/logs/${id}.out -eo ${runDir}/logs/${id}.err -W 12:00 -n 1 -R "rusage[mem=8000]" sh ${runDir}/logs/${id}-lumpy.sh

       #rm -rf ${id}.discordants.unsorted.bam ${id}.splitters.unsorted.bam ${id}.discordants.bam ${id}.splitters.bam

    fi

done < $bamFileList

