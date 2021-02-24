#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_Manta.sh
## Usage: ./collect_02_Manta.sh myConfigFile.config

## This script must be ran after all jobs from run_02_Manta.sh are finished.
## A folder Results will be created inside the folder 02_SVtools/Manta.

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
##
# MANTA
# DUP, DEL, INS, INV, BND
# Filters:: 
# 01: Keep only calls from chr 1-22+X+Y
# 02: Removing SVs smaller than 50 bp
# 03: Keep only PASS calls
##

## Run folder
mkdir -p ${outputFolderPath} || exit 1
prevResultsDir=${outputFolderPath}"/02_SVtools/Manta"
runDir=${prevResultsDir}"/Results"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/svCalls || exit 1
mkdir -p ${runDir}/svFiltered || exit 1
mkdir -p ${runDir}/svBenchmark || exit 1
mkdir -p ${runDir}/SURVIVOR || exit 1

touch ${outputFolderPath}/svCalls.blacklist

## Loading modules
module purge
module load R
module load vcftools
module load bcftools

# =============================================================================
#                           Collect results
# =============================================================================

## First check if all samples were processed and resulting files exist
echo "List with bam files: $bamFileList"
total_num_samples=$(cat $bamFileList | wc -l)
echo "Number of samples in the list: $total_num_samples"
COUNTER=0
while [ 1 ]
do

    read bamFile || break
    fileName=$(basename "$bamFile")
    filename=$(basename -- "$bamFile")
    extension="${filename##*.}"
    filename="${filename%.*}"
    id=`echo $fileName | cut -d '.' -f 1`

    if [ ! -s "${prevResultsDir}/data/${id}/results/variants/diploidSV.vcf.gz" ]; then
       COUNTER=$[$COUNTER +1]
       # Add to blacklist
       print "${id}\tMANTA_failed" >> ${outputFolderPath}/svCalls.blacklist
    fi

done < $bamFileList
echo "Number of failed samples: $COUNTER"
totalFails=$((COUNTER))
if [ $totalFails -gt 0 ]; then
    echo "Some samples had failed. Please fix first"
    exit 1
fi

## Loop over each sample
while [ 1 ]
do

    read bamFile || break

    fileName=$(basename "$bamFile")
    filename=$(basename -- "$bamFile")
    extension="${filename##*.}"

    fileName=$(basename "$bamFile")
    id=`echo $fileName | cut -d '.' -f 1`
    echo "|-- $bamFile ---- $id --|"
    
    if [ ! -s ${runDir}/svFiltered/${id}.vcf.filt.gz ]; then
    
       ## Copy files into the Results/svCalls folder
       cp ${prevResultsDir}/data/${id}/results/variants/diploidSV.vcf.gz  ${runDir}/svCalls/${id}.vcf.gz
       cp ${prevResultsDir}/data/${id}/results/variants/diploidSV.vcf.gz.tbi  ${runDir}/svCalls/${id}.vcf.gz.tbi
       # cp ${prevResultsDir}/data/${id}/results/variants/candidateSV.vcf.gz  ${runDir}/${id}-candidateSV.vcf.gz
       # cp ${prevResultsDir}/data/${id}/results/variants/candidateSV.vcf.gz.tbi  ${runDir}/${id}-candidateSV.vcf.gz.tbi
       # cp ${prevResultsDir}/data/${id}/results/variants/candidateSmallIndels.vcf.gz  ${runDir}/${id}-candidateSmallIndels.vcf.gz
       # cp ${prevResultsDir}/data/${id}/results/variants/candidateSmallIndels.vcf.gz.tbi  ${runDir}/${id}-candidateSmallIndels.vcf.gz.tbi

       ## Files are originally compressed. So we first decompress them and sort
       zcat ${runDir}/svCalls/${id}.vcf.gz > ${runDir}/svCalls/${id}.vcf
       vcf-sort -c ${runDir}/svCalls/${id}.vcf > ${runDir}/svCalls/${id}.sorted.vcf
       mv ${runDir}/svCalls/${id}.sorted.vcf ${runDir}/svCalls/${id}.vcf

       ## Clean and organize VCF
       # Filter 01: keep only calls from chr 1-22+X+Y and Remove partially-assembled (ie. inferred) insertions
       bcftools view ${runDir}/svCalls/${id}.vcf.gz --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y | grep -v "<INS>" > ${runDir}/svFiltered/${id}.vcf.01

       # Filter 02: Removing SVs smaller than 50 bp.
       SURVIVOR filter ${runDir}/svFiltered/${id}.vcf.01 NA 50 -1 0 -1 ${runDir}/svFiltered/${id}.vcf.02

       # Filter 03: Keep only PASS calls.
       awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${runDir}/svFiltered/${id}.vcf.02 > ${runDir}/svFiltered/${id}.vcf.03

       # Rename sample id
       ${svpipeline_r_lib}/rename_SampleNamesInVCF.R ${runDir}/svFiltered/${id}.vcf.03 MANTA_${id} ${runDir}/svFiltered/${id}.vcf.04

       # Rename REF and ALF fields
       ${svpipeline_r_lib}/rename_REFandALT.R ${runDir}/svFiltered/${id}.vcf.04 ${runDir}/svFiltered/${id}.vcf.05

       # Rename SV ids
       ${svpipeline_r_lib}/rename_SVIds.R ${runDir}/svFiltered/${id}.vcf.05 ${runDir}/svFiltered/${id}.vcf.06.gz
       gunzip -f ${runDir}/svFiltered/${id}.vcf.06.gz # decompress

       # Copy/Rename last step to .filt
       cp ${runDir}/svFiltered/${id}.vcf.06 ${runDir}/svFiltered/${id}.vcf.filt

       # Extract Genotypes for each SV
       bgzip -c ${runDir}/svFiltered/${id}.vcf.filt > ${runDir}/svFiltered/${id}.vcf.filt.gz
       tabix -p vcf ${runDir}/svFiltered/${id}.vcf.filt.gz
       bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${runDir}/svFiltered/${id}.vcf.filt.gz > ${runDir}/svFiltered/${id}.vcf.filt.genotypes
    
       # For Truvari compatibility - If your vcf has DUP calls that generally indicate insertions of tandem duplications (e.g., manta), then convert SVTYPE=DUP to SVTYPE=INS
       sed 's/SVTYPE=DUP/SVTYPE=INS/' ${runDir}/svFiltered/${id}.vcf.filt | grep -E '^#|SVTYPE=INS' > ${runDir}/svFiltered/${id}.vcf.filt.DUPtoINS

       # Splitting VCF into DEL, INS, DUP, INV, BND
       grep -E '^#|SVTYPE=DEL' ${runDir}/svFiltered/${id}.vcf.filt > ${runDir}/svFiltered/${id}.vcf.filt.DEL

       grep -E '^#|SVTYPE=DUP' ${runDir}/svFiltered/${id}.vcf.filt > ${runDir}/svFiltered/${id}.vcf.filt.DUP

       grep -E '^#|SVTYPE=INS' ${runDir}/svFiltered/${id}.vcf.filt > ${runDir}/svFiltered/${id}.vcf.filt.INS

       grep -E '^#|SVTYPE=INV' ${runDir}/svFiltered/${id}.vcf.filt > ${runDir}/svFiltered/${id}.vcf.filt.INV

       grep -E '^#|SVTYPE=BND' ${runDir}/svFiltered/${id}.vcf.filt > ${runDir}/svFiltered/${id}.vcf.filt.BND
    
       # BNDs are actually TRAs? 
       sed 's/SVTYPE=BND/SVTYPE=TRA/' ${runDir}/svFiltered/${id}.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${runDir}/svFiltered/${id}.vcf.filt.BNDtoTRA

    fi

done < $bamFileList

## SUBMIT TO BENCHMARK (only for DEL, INS and DUPtoINS)
while [ 1 ]
do

    read bamFile || break

    fileName=$(basename "$bamFile")
    filename=$(basename -- "$bamFile")
    extension="${filename##*.}"

    fileName=$(basename "$bamFile")
    id=`echo $fileName | cut -d '.' -f 1`
    echo "|-- $bamFile ---- $id --|"

    # Prepare script for benchmarking
    cat /dev/null > "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#!/bin/bash" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -L /bin/bash" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -n 1" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -R rusage[mem=4000]" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -J ${id}.bench" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -oo /dev/null" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -eo /dev/null" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -q premium" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -W 144:00" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "#BSUB -P acc_ad-omics" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "${bench_truvari} ${runDir}/svCalls/${id}.vcf ${runDir}/svBenchmark/${id}_raw" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt ${runDir}/svBenchmark/${id}_filt" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt.DUPtoINS ${runDir}/svBenchmark/${id}_filt.DUPtoINS" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt.DEL ${runDir}/svBenchmark/${id}_filt.DEL" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    echo "${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt.INS ${runDir}/svBenchmark/${id}_filt.INS" >> "${runDir}/svBenchmark/${id}_benchmark.sh"
    chmod a+x "${runDir}/svBenchmark/${id}_benchmark.sh"
    if [ ! -s ${runDir}/svBenchmark/${id}_filt/${id}.giab_report.txt ]; then
       echo "Report for ${id}"
       #bsub < ${runDir}/svBenchmark/${id}_benchmark.sh
       ${bench_truvari} ${runDir}/svCalls/${id}.vcf ${runDir}/svBenchmark/${id}_raw
       ${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt ${runDir}/svBenchmark/${id}_filt
       ${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt.DUPtoINS ${runDir}/svBenchmark/${id}_filt.DUPtoINS
       ${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt.DEL ${runDir}/svBenchmark/${id}_filt.DEL
       ${bench_truvari} ${runDir}/svFiltered/${id}.vcf.filt.INS ${runDir}/svBenchmark/${id}_filt.INS
    fi

done < $bamFileList

# =============================================================================
#                           Merge population call set
# =============================================================================

# Merging samples with SURVIVOR
breakpoint_dist=1000 # max distance between breakpoints
min_num_calls=1 # Minimum number of supporting sample
use_type=1 # Take the type into account (1==yes, else no)
use_strand=1 # Take the strands of SVs into account (1==yes, else no)
dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
min_sv_size=50 # Minimum size of SVs to be taken into account.

find ${runDir}/svFiltered -name "*.vcf.filt" -type f > ${runDir}/SURVIVOR/vcfCallFiles_ALL.list
SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_ALL.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_ALL.vcf
vcf-sort -c ${runDir}/SURVIVOR/callsMerged_ALL.vcf > ${runDir}/SURVIVOR/callsMerged_ALL.sorted.vcf
bgzip -c ${runDir}/SURVIVOR/callsMerged_ALL.sorted.vcf > ${runDir}/SURVIVOR/callsMerged_ALL.sorted.vcf.gz
tabix -p vcf ${runDir}/SURVIVOR/callsMerged_ALL.sorted.vcf.gz
bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO ${runDir}/SURVIVOR/callsMerged_ALL.sorted.vcf.gz > ${runDir}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf
bgzip -c ${runDir}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf > ${runDir}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf.gz
tabix -p vcf ${runDir}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf.gz

find ${runDir}/svFiltered -name "*.vcf.filt.INS" -type f > ${runDir}/SURVIVOR/vcfCallFiles_INS.list
SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_INS.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_INS.vcf

find ${runDir}/svFiltered -name "*.vcf.filt.DEL" -type f > ${runDir}/SURVIVOR/vcfCallFiles_DEL.list
SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_DEL.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_DEL.vcf

find ${runDir}/svFiltered -name "*.vcf.filt.DUP" -type f > ${runDir}/SURVIVOR/vcfCallFiles_DUP.list
SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_DUP.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_DUP.vcf

find ${runDir}/svFiltered -name "*.vcf.filt.INV" -type f > ${runDir}/SURVIVOR/vcfCallFiles_INV.list
SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_INV.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_INV.vcf

find ${runDir}/svFiltered -name "*.vcf.filt.BND" -type f > ${runDir}/SURVIVOR/vcfCallFiles_BND.list
SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_BND.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_BND.vcf

