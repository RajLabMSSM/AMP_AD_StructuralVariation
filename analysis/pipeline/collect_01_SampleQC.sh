#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 01_SampleQC.sh
## Usage: ./collect_01_SampleQC.sh myConfigFile.config

## This script must be ran after all jobs from run_01_SampleQC.sh are finished.
## A folder Results will be created inside the folder 01_SampleQC.
## Inside this folder a file (cohort.QC.metrics) containing a summary QC metrics of all samples will be created.

## This step is based on the the first module from the HOLMES pipeline. Some of its scripts are included in the QC_Scripts folder.
## Collins RL, et al. Defining the spectrum of large inversions, complex structural variation, and chromothripsis in the morbid genome. Genome Biol. (2017)

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
prevResultsDir=${outputFolderPath}"/01_SampleQC"
runDir=${outputFolderPath}"/01_SampleQC/Results"
mkdir -p ${runDir} || exit 1

# Loading modules
module load sambamba
module load picard
module load samtools
module load bamtools

# =============================================================================
#                   Check if all files QC metrics were collected
# =============================================================================

while [ 1 ]
do
    read line || break
    name=$(basename "$line")
    id=`echo $name | cut -d '.' -f 1`
    if [ -z "${prevResultsDir}/${id}.aneuploidyCheck" ]; then
      echo "Metric failed for ${id}.aneuploidyCheck. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.flagstat" ]; then
      echo "Metric failed for ${id}.flagstat. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.metrics.alignment_summary_metrics" ]; then
      echo "Metric failed for ${id}.metrics.alignment_summary_metrics. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.metrics.insert_size_metrics" ]; then
      echo "Metric failed for ${id}.metrics.insert_size_metrics. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.sexCheck" ]; then
      echo "Metric failed for ${id}.sexCheck. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.stats" ]; then
      echo "Metric failed for ${id}.stats. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.wgs" ]; then
      echo "Metric failed for ${id}.wgs. Please check before running again."
      exit 1
    fi
    if [ -z "${prevResultsDir}/${id}.complexity" ]; then
      echo "Metric failed for ${id}.complexity. Please check before running again."
      exit 1
    fi
done < $bamFileList

# =============================================================================
#                            Gather results
# =============================================================================

## Parse aneuploidy check
paste <( echo -e "chr\texpected" ) <( cut -f1 ${bamFileList} | sed 's!.*/!!' | cut -d '.' -f 1 | paste -s ) > ${runDir}/cohort.aneuploidyCheck.fractions.txt
paste <( echo -e "chr\texpected" ) <( cut -f1 ${bamFileList} | sed 's!.*/!!' | cut -d '.' -f 1 |paste -s ) > ${runDir}/cohort.aneuploidyCheck.copies.txt
echo -e "$( seq 1 22 )\nX\nY" > ${runDir}/build.fractions.tmp
echo -e "$( seq 1 22 )\nX\nY" > ${runDir}/build.copies.tmp
ID=$( head -n1 ${bamFileList} | cut -f1 )
name=$(basename $ID)
id=`echo $name | cut -d '.' -f 1`
paste ${runDir}/build.fractions.tmp <( fgrep -v "#" ${prevResultsDir}/${id}.aneuploidyCheck | sed '/^$/d' | cut -f5 ) > ${runDir}/build.fractions.tmp2
paste ${runDir}/build.copies.tmp <( perl -E "say \"2\n\" x 24" | sed '/^$/d' ) > ${runDir}/build.copies.tmp2
mv ${runDir}/build.fractions.tmp2 ${runDir}/build.fractions.tmp
mv ${runDir}/build.copies.tmp2 ${runDir}/build.copies.tmp
while [ 1 ]
do
  read line || break
  name=$(basename "$line")
  id=`echo $name | cut -d '.' -f 1`
  paste ${runDir}/build.fractions.tmp <( fgrep -v "#" ${prevResultsDir}/${id}.aneuploidyCheck | sed '/^$/d' | cut -f4 ) > ${runDir}/build.fractions.tmp2
  paste ${runDir}/build.copies.tmp <( fgrep -v "#" ${prevResultsDir}/${id}.aneuploidyCheck | sed '/^$/d' | cut -f6 ) > ${runDir}/build.copies.tmp2
  mv ${runDir}/build.fractions.tmp2 ${runDir}/build.fractions.tmp
  mv ${runDir}/build.copies.tmp2 ${runDir}/build.copies.tmp
done < ${bamFileList}
cat ${runDir}/build.fractions.tmp >> ${runDir}/cohort.aneuploidyCheck.fractions.txt
cat ${runDir}/build.copies.tmp >> ${runDir}/cohort.aneuploidyCheck.copies.txt
rm ${runDir}/build.fractions.tmp ${runDir}/build.copies.tmp

## Collect summary for each sample
echo -e "ID\tTotal_Reads\tTotal_Reads_Mapped\tTotal_Reads_Mapped_Perc\tReads_Fwd\tReads_Rev\tFailed_QC\tReads_Duplicated\tPaired_end\tProper_pairs\tBoth_pairs\tSingletons\tAvg_Insert_Size\tMedian_Insert_Size\tPair_Dup\tHap_Nuc_Cov\tPct_Chimera\tAvg_Read_Length\tPct_Adapter\tObserved_Sex" > ${runDir}/cohort.QC.metrics
while [ 1 ]
do
  read line || break
  name=$(basename "$line")
  id=`echo $name | cut -d '.' -f 1`
  total=$( grep '^Total reads:' ${prevResultsDir}/${id}.stats | awk '{ print $3 }' ) # total reads in BAM
  rd_mapped=$( grep '^Mapped reads:' ${prevResultsDir}/${id}.stats | awk '{ print $3 }' ) # reads mapped
  rd_mapped_perc=$( grep '^Mapped reads:' ${prevResultsDir}/${id}.stats | awk '{ print $4 }' | sed -e 's/(//g ; s/)//g ; s/%//g' ) # read mapped %
  rd_fwr=$( grep '^Forward strand:' ${prevResultsDir}/${id}.stats | awk '{ print $3 }' ) # reads forward strand mapped
  rd_rev=$( grep '^Reverse strand:' ${prevResultsDir}/${id}.stats | awk '{ print $3 }' ) # reads reverse strand mapped
  failed_qc=$( grep '^Failed QC:' ${prevResultsDir}/${id}.stats | awk '{ print $3 }' ) # Failed QC
  rd_dup=$( grep '^Duplicates:' ${prevResultsDir}/${id}.stats | awk '{ print $2 }' ) # reads duplicated
  rd_pe=$( grep '^Paired-end reads:' ${prevResultsDir}/${id}.stats | awk '{ print $3 }' ) # paired-end reads
  prop=$( grep 'Proper-pairs' ${prevResultsDir}/${id}.stats | awk '{ print $2 }' ) # 'Proper-pairs'
  both_pe=$( grep '^Both pairs mapped:' ${prevResultsDir}/${id}.stats | awk '{ print $4 }' ) # Both pairs mapped
  singletons=$( grep '^Singletons:' ${prevResultsDir}/${id}.stats | awk '{ print $2 }' ) # Singletons
  avg_ins=$( grep '^Average insert size (absolute value):' ${prevResultsDir}/${id}.stats | awk '{ print $6 }' ) # Average insert size (absolute value)
  med_ins=$( grep '^Median insert size (absolute value):' ${prevResultsDir}/${id}.stats | awk '{ print $6 }' ) # Median insert size (absolute value)
  pr_dup=$( grep -A1 '^LIBRARY' ${prevResultsDir}/${id}.complexity | tail -n1 | awk '{ print $(NF-1) }' ) # pair dup rate
  ncov=$( grep -A1 '^GENOME_TERRITORY' ${prevResultsDir}/${id}.wgs | tail -n1 | awk '{ print $2 }' ) # mean nucleotide cov
  chim_rt=$( grep -A2 '^CATEGORY' ${prevResultsDir}/${id}.metrics.alignment_summary_metrics | tail -n2 | awk '{ print $21 }' | awk '{sum+=$0} END { print sum/NR}') # mean pct_chimeras (first of pairs and second of pair)
  read_len=$( grep -A2 '^CATEGORY' ${prevResultsDir}/${id}.metrics.alignment_summary_metrics | tail -n2 | awk '{ print $16 }' | awk '{sum+=$0} END { print sum/NR}') # mean read length (first of pairs and second of pair)
  adapter_rt=$( grep -A2 '^CATEGORY' ${prevResultsDir}/${id}.metrics.alignment_summary_metrics | tail -n2 | awk '{ print $22 }' | awk '{sum+=$0} END { print sum/NR}') # mean adapter pct (first of pairs and second of pair)
  osex=$( fgrep -w "PREDICTED SEX" ${prevResultsDir}/${id}.sexCheck | awk '{ print $3 }' ) # predicted sex from sexcheck.sh
  echo -e "${id}\t${total}\t${rd_mapped}\t${rd_mapped_perc}\t${rd_fwr}\t${rd_rev}\t${failed_qc}\t${rd_dup}\t${rd_pe}\t${prop}\t${both_pe}\t${singletons}\t${avg_ins}\t${med_ins}\t${pr_dup}\t${ncov}\t${chim_rt}\t${read_len}\t${adapter_rt}\t${osex}" #print metrics
done < ${bamFileList} >> ${runDir}/cohort.QC.metrics

