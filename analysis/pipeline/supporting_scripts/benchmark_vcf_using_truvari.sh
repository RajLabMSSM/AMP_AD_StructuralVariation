#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

module load vcftools
module load bcftools
module load samtools

if [ -d "$1" ]; then
   input=$(ls $1/*.vcf) # a folder with vcf files as input
else
   input=$1 # .vcf file as input
fi

output=$2

prepare_vcf_files()
{
  dir_path=$1
  for filename in ${dir_path}/*; do
    name=$(basename "$filename")
    id=`echo $name | cut -d '.' -f 1`
    echo "Sample: ${id}"
    echo "Sorting..."
    vcf-sort -c $filename > ${dir_path}/${id}.sorted.vcf 2> /dev/null
    echo "Compressing (bgzip)"
    bgzip -c ${dir_path}/${id}.sorted.vcf > ${dir_path}/${id}.sorted.vcf.gz
    echo "Indexing (tabix)"
    tabix -p vcf ${dir_path}/${id}.sorted.vcf.gz
    echo "Collecting status with SURVIVOR"
    SURVIVOR stats ${dir_path}/${id}.sorted.vcf -1 -1 -1 ${dir_path}/${id}.survivor.stats &> ${dir_path}/${id}.survivor.stats.summary
    echo "Benchmark (truvari)"
    rm -rf ${dir_path}/truvari/
    truvari \
       -b ${benchmark_vcf} \
       -c ${dir_path}/${id}.sorted.vcf.gz \
       -o ${dir_path}/truvari/ \
       -f ${referenceFasta} \
       --includebed ${benchmark_bed} \
       --passonly \
       --pctsim=0 \
       -r 2000 --giabreport
    rm -rf ${dir_path}/truvari_DEL/
    truvari \
       -b ${benchmark_vcf_DEL} \
       -c ${dir_path}/${id}.sorted.vcf.gz \
       -o ${dir_path}/truvari_DEL/ \
       -f ${referenceFasta} \
       --includebed ${benchmark_bed} \
       --passonly \
       --pctsim=0 \
       -r 2000 --giabreport
    rm -rf ${dir_path}/truvari_INS/
    truvari \
       -b ${benchmark_vcf_INS} \
       -c ${dir_path}/${id}.sorted.vcf.gz \
       -o ${dir_path}/truvari_INS/ \
       -f ${referenceFasta} \
       --includebed ${benchmark_bed} \
       --passonly \
       --pctsim=0 \
       -r 2000 --giabreport
    cp ${dir_path}/${id}.survivor.stats ${dir_path}/${id}.survivor.stats.summary ${dir_path}/${id}.survivor.stats_CHR ${dir_path}/../
    cp ${dir_path}/truvari/summary.txt ${dir_path}/../${id}.truvari
    cp ${dir_path}/truvari_DEL/summary.txt ${dir_path}/../${id}.truvari_DEL
    cp ${dir_path}/truvari_INS/summary.txt ${dir_path}/../${id}.truvari_INS
    cp ${dir_path}/truvari/giab_report.txt ${dir_path}/../${id}.giab_report.txt
    cp ${dir_path}/truvari_DEL/giab_report.txt ${dir_path}/../${id}.giab_report_DEL
    cp ${dir_path}/truvari_INS/giab_report.txt ${dir_path}/../${id}.giab_report_INS
  done
}

if [ -d "$output" ]; then
  # Control will enter here if $output exists.
  echo "Output directory already exists. Cleaning first."
  rm -rf $output
  # exit 1
fi

mkdir -p $output/ || exit 1
mkdir -p $output/temp || exit 1
cp $input $output/temp
prepare_vcf_files "$output/temp"

echo "Cleaning files."
rm -rf $output/temp $output/truvari*

echo "Done."
