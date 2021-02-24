#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 02_SVE.sh
## Usage: ./collect_02_SVE.sh myConfigFile.config

## This script must be ran after all jobs from run_02_SVE.sh are finished.
## A folder Results will be created inside the folder 02_SVtools/SVEPipeline.

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
# SVE pipeline tools
# 
# Filters:: 
# 01: Keep only calls from chr 1-22+X+Y
# 02: Removing SVs smaller than 50 bp
# 03: Keep only PASS calls
##

## Run folder
mkdir -p ${outputFolderPath} || exit 1
prevResultsDir=${outputFolderPath}"/02_SVtools/SVEPipeline"
runDir=${prevResultsDir}"/Results"
mkdir -p ${runDir} || exit 1

touch ${outputFolderPath}/svCalls.blacklist

## Loading modules
module purge
module load vcftools

# =============================================================================
#                           Collect results
# =============================================================================

## First check if all samples were processed and resulting files exist
echo "List with bam files: $bamFileList"
total_num_samples=$(cat $bamFileList | wc -l)
echo "Number of samples in the list: $total_num_samples"
BREAKDANCERCOUNTER=0 #S_4
BREAKSEQCOUNTER=0 #S_35
CNVNATORCOUNTER=0 #S_10
while [ 1 ]
do

    read bamFile || break
    fileName=$(basename "$bamFile")
    filename=$(basename -- "$bamFile")
    extension="${filename##*.}"
    filename="${filename%.*}" # File name with no extension
    id=`echo $fileName | cut -d '.' -f 1` # Clear file name, remove everything after '.'
    filename=$id

    if [ ! -s "${prevResultsDir}/svCalls/${filename}_S4.vcf" ]; then
       BREAKDANCERCOUNTER=$[$BREAKDANCERCOUNTER +1]
       # Add to blacklist
       printf "${id}\tBREAKDANCER_failed\n" >> ${outputFolderPath}/svCalls.blacklist
    fi

    if [ ! -s "${prevResultsDir}/svCalls/${filename}_S35.vcf" ]; then
       BREAKSEQCOUNTER=$[$BREAKSEQCOUNTER +1]
       # Add to blacklist
       printf "${id}\tBREAKSEQ_failed\n" >> ${outputFolderPath}/svCalls.blacklist
    fi

    if [ ! -s "${prevResultsDir}/svCalls/${filename}_S10.vcf" ]; then
       CNVNATORCOUNTER=$[$CNVNATORCOUNTER +1]
       # Add to blacklist
       printf "${id}\tCNVNATOR_failed\n" >> ${outputFolderPath}/svCalls.blacklist
    fi

done < $bamFileList

echo "Number of failed samples (BreakDancer): $BREAKDANCERCOUNTER"
echo "Number of failed samples (BreakSeq): $BREAKSEQCOUNTER"
echo "Number of failed samples (CNVnator): $CNVNATORCOUNTER"
totalFails=$((BREAKDANCERCOUNTER + BREAKSEQCOUNTER + CNVNATORCOUNTER))
if [ $totalFails -gt 0 ]; then
    echo "Some samples had failed. Please fix first"
    # exit 1
fi

run_collect() 
{
  sve_id=$1
  tool_name=$2
  while [ 1 ]
  do

      read bamFile || break

      filename=$(basename -- "$bamFile")
      extension="${filename##*.}"
      filename="${filename%.*}"

      fileName=$(basename "$bamFile")
      id=`echo $fileName | cut -d '.' -f 1`

      echo "|-- $bamFile ---- $id --|"
  
      if [ -n "${prevResultsDir}/svCalls/${filename}_${sve_id}" ]; then
        
        if [[ -s "${prevResultsDir}/svCalls/${filename}_${sve_id}" ]]; then # Workaround for filenames differences (id.final.bam vs id.bam). Still need to be fixed properly. 
            cp ${prevResultsDir}/svCalls/${filename}_${sve_id}  ${runDir}/${tool_name}/svCalls/${id}.vcf
        else
            cp ${prevResultsDir}/svCalls/${id}_${sve_id}  ${runDir}/${tool_name}/svCalls/${id}.vcf
        fi

        if [ ! -s ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt ]; then

           rm -rf ${runDir}/${tool_name}/svFiltered/${id}.*          
    
           # Sort vcf
           vcf-sort -c ${runDir}/${tool_name}/svCalls/${id}.vcf > ${runDir}/${tool_name}/svCalls/${id}.sorted.vcf
           mv ${runDir}/${tool_name}/svCalls/${id}.sorted.vcf ${runDir}/${tool_name}/svCalls/${id}.vcf

           # Compress and index vcf
           bgzip -c ${runDir}/${tool_name}/svCalls/${id}.vcf > ${runDir}/${tool_name}/svCalls/${id}.vcf.gz
           tabix -p vcf ${runDir}/${tool_name}/svCalls/${id}.vcf.gz

           ## Clean and organize VCF
           # Filter 01: keep only calls from chr 1-22+X+Y
           bcftools view ${runDir}/${tool_name}/svCalls/${id}.vcf.gz --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y > ${runDir}/${tool_name}/svFiltered/${id}.vcf.01

           # Filter 02: Removing SVs smaller than 50 bp.
           SURVIVOR filter ${runDir}/${tool_name}/svFiltered/${id}.vcf.01 NA 50 -1 10000000 -1 ${runDir}/${tool_name}/svFiltered/${id}.vcf.02

           # Filter 03: Keep only PASS calls - already filtered!
           # awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${runDir}/svFiltered/${id}.vcf.02 > ${runDir}/svFiltered/${id}.vcf.03

           # Rename sample id
           ${svpipeline_r_lib}/rename_SampleNamesInVCF.R ${runDir}/${tool_name}/svFiltered/${id}.vcf.02 ${tool_name}_${id} ${runDir}/${tool_name}/svFiltered/${id}.vcf.03

           # Rename REF and ALF fields
           ${svpipeline_r_lib}/rename_REFandALT.R ${runDir}/${tool_name}/svFiltered/${id}.vcf.03 ${runDir}/${tool_name}/svFiltered/${id}.vcf.04.gz
           gunzip ${runDir}/${tool_name}/svFiltered/${id}.vcf.04.gz # decompress

           # Copy/Rename last step to .filt
           cp ${runDir}/${tool_name}/svFiltered/${id}.vcf.04 ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt

           # Extract Genotypes for each SV
           bgzip -c ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.gz
           tabix -p vcf ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.gz
           bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.gz > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.genotypes

           # For Truvari compatibility - If your vcf has DUP calls that generally indicate insertions of tandem duplications (e.g., manta), then convert SVTYPE=DUP to SVTYPE=INS
           sed 's/SVTYPE=DUP/SVTYPE=INS/' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt | grep -E '^#|SVTYPE=INS' > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.DUPtoINS

           # Splitting VCF into DEL, INS, DUP, INV, TRA, BND
           grep -E '^#|SVTYPE=DEL' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.DEL

           grep -E '^#|SVTYPE=DUP' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.DUP

           grep -E '^#|SVTYPE=INS' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.INS

           grep -E '^#|SVTYPE=INV' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.INV

           grep -E '^#|SVTYPE=BND' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.BND

           grep -E '^#|SVTYPE=TRA' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.TRA

           # BNDs are actually TRAs? 
           sed 's/SVTYPE=BND/SVTYPE=TRA/' ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.BNDtoTRA

        fi
      fi

      # Prepare script for benchmarking
      cat /dev/null > "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      echo "${bench_truvari} ${runDir}/${tool_name}/svCalls/${id}.vcf ${runDir}/${tool_name}/svBenchmark/${id}_raw" >> "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      echo "${bench_truvari} ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt ${runDir}/${tool_name}/svBenchmark/${id}_filt" >> "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      echo "${bench_truvari} ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.DUPtoINS ${runDir}/${tool_name}/svBenchmark/${id}_filt.DUPtoINS" >> "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      echo "${bench_truvari} ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.DEL ${runDir}/${tool_name}/svBenchmark/${id}_filt.DEL" >> "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      echo "${bench_truvari} ${runDir}/${tool_name}/svFiltered/${id}.vcf.filt.INS ${runDir}/${tool_name}/svBenchmark/${id}_filt.INS" >> "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      chmod a+x "${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh"
      bsub -n 1 -L /bin/bash -R "rusage[mem=6000]" -W 144:00 -oo /dev/null -eo /dev/null -m bode -P acc_ad-omics -q premium -J ${id}.bench ${runDir}/${tool_name}/svBenchmark/${id}_benchmark.sh

  done < ${bamFileList}

  # Merging samples with SURVIVOR
  breakpoint_dist=1000 # max distance between breakpoints
  min_num_calls=1 # Minimum number of supporting sample
  use_type=1 # Take the type into account (1==yes, else no)
  use_strand=1 # Take the strands of SVs into account (1==yes, else no)
  dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
  min_sv_size=50 # Minimum size of SVs to be taken into account.

  find ${runDir}/${tool_name}/svFiltered -name "*.vcf.filt" -type f > ${runDir}/${tool_name}/SURVIVOR/vcfCallFiles_ALL.list
  SURVIVOR merge ${runDir}/${tool_name}/SURVIVOR/vcfCallFiles_ALL.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.vcf
  vcf-sort -t ${tmpDir} -c ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.vcf > ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.vcf
  bgzip -c ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.vcf > ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.vcf.gz
  tabix -p vcf ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.vcf.gz
  bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.vcf.gz > ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf
  bgzip -c ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf > ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf.gz
  tabix -p vcf ${runDir}/${tool_name}/SURVIVOR/callsMerged_ALL.sorted.clean.vcf.gz

  find ${runDir}/${tool_name}/svFiltered -name "*.vcf.filt.INS" -type f > ${runDir}/${tool_name}/SURVIVOR/vcfCallFiles_INS.list
  SURVIVOR merge ${runDir}/${tool_name}/SURVIVOR/vcfCallFiles_INS.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/${tool_name}/SURVIVOR/callsMerged_INS.vcf

  find ${runDir}/${tool_name}/svFiltered -name "*.vcf.filt.DEL" -type f > ${runDir}/${tool_name}/SURVIVOR/vcfCallFiles_DEL.list
  SURVIVOR merge ${runDir}/${tool_name}/SURVIVOR/vcfCallFiles_DEL.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/${tool_name}/SURVIVOR/callsMerged_DEL.vcf

  find ${runDir}/svFiltered -name "*.vcf.filt.DUP" -type f > ${runDir}/SURVIVOR/vcfCallFiles_DUP.list
  SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_DUP.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_DUP.vcf

  find ${runDir}/svFiltered -name "*.vcf.filt.INV" -type f > ${runDir}/SURVIVOR/vcfCallFiles_INV.list
  SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_INV.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_INV.vcf

  find ${runDir}/svFiltered -name "*.vcf.filt.BND" -type f > ${runDir}/SURVIVOR/vcfCallFiles_BND.list
  SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_BND.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_BND.vcf

  find ${runDir}/svFiltered -name "*.vcf.filt.TRA" -type f > ${runDir}/SURVIVOR/vcfCallFiles_TRA.list
  SURVIVOR merge ${runDir}/SURVIVOR/vcfCallFiles_TRA.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${runDir}/SURVIVOR/callsMerged_TRA.vcf

}

echo "######################################################"
echo "#########  Gathering results for each tool..."

## BreakDancer
## SVE code: S4
## DEL, DUP, INV, TRA
if ls ${prevResultsDir}/svCalls/*_S4.vcf>/dev/null 2>&1; then
   ls ${prevResultsDir}/svCalls/*_S4.vcf > ${runDir}/BreakDancer_S4_samples_vcf.list
   numfiles=$(wc -l ${runDir}/BreakDancer_S4_samples_vcf.list | awk '{print $1}')
   echo "BreakDancer: ${numfiles} samples (.vcf) found."
   mkdir -p ${runDir}/BreakDancer || exit 1
   mkdir -p ${runDir}/BreakDancer/svCalls || exit 1
   mkdir -p ${runDir}/BreakDancer/svFiltered || exit 1
   mkdir -p ${runDir}/BreakDancer/svBenchmark || exit 1
   mkdir -p ${runDir}/BreakDancer/SURVIVOR || exit 1
   # Adding FORMAT and Genotype columns for calls using only 1 sample 
   while [ 1 ]
   do
      read vcfFile || break
      name=$(basename "$vcfFile")
      id=`echo $name | cut -d '.' -f 1`
      cp ${vcfFile} ${vcfFile}.original
      cat /dev/null > ${vcfFile}.fixColumn
      grep "^##" ${vcfFile} >> ${vcfFile}.fixColumn
      echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${vcfFile}.fixColumn
      echo '##FORMAT=<ID=CN,Number=1,Type=String,Description="Copy number genotype for imprecise events">' >> ${vcfFile}.fixColumn
      # Add FORMAT and Sample to Header
      grep "^#CHROM" ${vcfFile} | awk -v newdata="FORMAT\t${id}" 'BEGIN{FS=OFS="\t"} {print $0 OFS newdata}' >> ${vcfFile}.fixColumn
      # Add format and genotype to calls (generic)
      grep -v "^#" ${vcfFile} | awk -v newdata="GT:CN\t./.:." 'BEGIN{FS=OFS="\t"} {print $0 OFS newdata}' >> ${vcfFile}.fixColumn
   done < ${runDir}/BreakDancer_S4_samples_vcf.list
   run_collect "S4.vcf.fixColumn" "BreakDancer"
else
   echo "No results for BreakDancer"
fi

## CNVnator
## SVE code: S10
## DEL, DUP
if ls ${prevResultsDir}/svCalls/*_S10.vcf >/dev/null 2>&1; then
   ls ${prevResultsDir}/svCalls/*_S10.vcf > ${runDir}/CNVnator_S10_samples_vcf.list
   echo "CNVnator: $(wc -l ${runDir}/CNVnator_S10_samples_vcf.list | awk '{print $1}') samples (.vcf) found."
   mkdir -p ${runDir}/CNVnator || exit 1
   mkdir -p ${runDir}/CNVnator/svCalls || exit 1
   mkdir -p ${runDir}/CNVnator/svFiltered || exit 1
   mkdir -p ${runDir}/CNVnator/svBenchmark || exit 1
   mkdir -p ${runDir}/CNVnator/SURVIVOR || exit 1
   # Keep PASS only
   while [ 1 ]
   do
      read vcfFile || break
      cp ${vcfFile} ${vcfFile}.original
      awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${vcfFile} > ${vcfFile}.pass
   done < ${runDir}/CNVnator_S10_samples_vcf.list
   run_collect "S10.vcf.pass" "CNVnator"
else
   echo "No results for CNVnator"
fi

## DELLY (use standalone instead)
## SV code: S11
if ls ${prevResultsDir}/svCalls/*_S11.vcf >/dev/null 2>&1; then
   ls ${prevResultsDir}/svCalls/*_S11.vcf > ${runDir}/DELLY_S11_samples_vcf.list
   echo "DELLY: $(wc -l ${runDir}/DELLY_S11_samples_vcf.list | awk '{print $1}') samples (.vcf) found."
   mkdir -p ${runDir}/DELLY || exit 1
   mkdir -p ${runDir}/DELLY/svCalls || exit 1
   mkdir -p ${runDir}/DELLY/svFiltered || exit 1
   mkdir -p ${runDir}/DELLY/svBenchmark || exit 1
   mkdir -p ${runDir}/DELLY/SURVIVOR || exit 1
   # Keep PASS only
   while [ 1 ]
   do
      read vcfFile || break
      cp ${vcfFile} ${vcfFile}.original
      awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${vcfFile} > ${vcfFile}.pass
   done < ${runDir}/DELLY_S11_samples_vcf.list
   run_collect "S11.vcf.pass" "DELLY"
else
   echo "No results for DELLY"
fi

## LUMPY (use standalone instead)
## SV code: S18
if ls ${prevResultsDir}/svCalls/*_S18.vcf >/dev/null 2>&1; then
   ls ${prevResultsDir}/svCalls/*_S18.vcf > ${runDir}/LUMPY_S18_samples_vcf.list
   echo "LUMPY: $(wc -l ${runDir}/LUMPY_S18_samples_vcf.list | awk '{print $1}') samples (.vcf) found."
   mkdir -p ${runDir}/LUMPY || exit 1
   mkdir -p ${runDir}/LUMPY/svCalls || exit 1
   mkdir -p ${runDir}/LUMPY/svFiltered || exit 1
   mkdir -p ${runDir}/LUMPY/svBenchmark || exit 1
   mkdir -p ${runDir}/LUMPY/SURVIVOR || exit 1
   # Keep PASS only
   while [ 1 ]
   do
      read vcfFile || break
      cp ${vcfFile} ${vcfFile}.original
      awk -F '\t' '{if($0 ~ /\#/) print; else {if($7 == ".") $7="PASS"; print} }' OFS='\t' ${vcfFile} > ${vcfFile}.pass
   done < ${runDir}/LUMPY_S18_samples_vcf.list
   run_collect "S18.vcf.pass" "LUMPY"
else
   echo "No results for LUMPY"
fi

## BreakSeq
## SV code: S35
## DEL, INS
if ls ${prevResultsDir}/svCalls/*_S35.vcf >/dev/null 2>&1; then
   ls ${prevResultsDir}/svCalls/*_S35.vcf > ${runDir}/BreakSeq_S35_samples_vcf.list
   echo "BreakSeq: $(wc -l ${runDir}/BreakSeq_S35_samples_vcf.list | awk '{print $1}') samples (.vcf) found."
   mkdir -p ${runDir}/BreakSeq || exit 1
   mkdir -p ${runDir}/BreakSeq/svCalls || exit 1
   mkdir -p ${runDir}/BreakSeq/svFiltered || exit 1
   mkdir -p ${runDir}/BreakSeq/svBenchmark || exit 1
   mkdir -p ${runDir}/BreakSeq/SURVIVOR || exit 1
   # Keep PASS only
   while [ 1 ]
   do
      read vcfFile || break
      cp ${vcfFile} ${vcfFile}.original
      awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' ${vcfFile} > ${vcfFile}.pass
   done < ${runDir}/BreakSeq_S35_samples_vcf.list
   run_collect "S35.vcf.pass" "BreakSeq"
else
   echo "No results for BreakSeq"
fi
