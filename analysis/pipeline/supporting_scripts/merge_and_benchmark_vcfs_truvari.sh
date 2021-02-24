#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

input=$1 # input is a file with a list of vcfs
output=$2 # output folder to write results
min_num_calls=$3 # min number of supporting caller
comp_against=$4 # options are "ALL", "DEL" or "INS"

module load vcftools
module load bcftools
module load samtools

# SURVIVOR inputs::
# File with VCF names and paths
breakpoint_dist=1000 # max distance between breakpoints
min_num_calls=1 # Minimum number of supporting caller
use_type=1 # Take the type into account (1==yes, else no)
use_strand=1 # Take the strands of SVs into account (1==yes, else no)
dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
min_sv_size=50 # Minimum size of SVs to be taken into account.

mkdir -p $output

# If input has only 1 vcf, merging can be skipped...
if [[ $(cat $input | sed '/^\s*$/d' | wc -l) -ge 1 ]]; then
   filename=${output}/merged.vcf
   # Merge vcfs using SURVIVOR
   SURVIVOR merge ${input} ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${filename}
   sed -i 's|FORMAT=<ID=DR,Number=1,Type=Integer|FORMAT=<ID=DR,Number=1,Type=String|g' ${filename}
   sed -i 's|ID=LN,Number=1,Type=Integer|ID=LN,Number=1,Type=String|g' ${filename}
   #sed -i 's|;AVGLEN=|;SVLEN=|g' ${filename}
else
   filename=$(cat ${input})
fi

# Extract Genotypes for each SV
echo "Sorting..."
vcf-sort -c ${filename} > ${filename}.sorted
mv ${filename}.sorted ${filename}
echo "Compressing (bgzip)"
bgzip -c ${filename} > ${filename}.gz
echo "Indexing (tabix)"
tabix -p vcf ${filename}.gz
echo "Saving genotypes"
bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${filename}.gz > ${filename}.genotypes

dir_path=${output}

name=$(basename "$filename")
id=`echo $name | cut -d '.' -f 1`
echo "Sample: ${id}"
#echo "Sorting..."
#vcf-sort -c $filename > ${dir_path}/${id}.sorted.vcf 2> /dev/null
#echo "Compressing (bgzip)"
#bgzip -c ${dir_path}/${id}.sorted.vcf > ${dir_path}/${id}.sorted.vcf.gz
#echo "Indexing (tabix)"
#tabix -p vcf ${dir_path}/${id}.sorted.vcf.gz
echo "Collecting status with SURVIVOR"
SURVIVOR stats ${dir_path}/${id}.vcf -1 -1 -1 ${dir_path}/${id}.survivor.stats &> ${dir_path}/${id}.survivor.stats.summary
echo "Benchmark (truvari)"
rm -rf ${dir_path}/truvari/
if [ "$comp_against" = "ALL" ]; then
   rm -rf ${dir_path}/truvari/
   truvari \
      -b ${benchmark_vcf} \
      -c ${dir_path}/${id}.vcf.gz \
      -o ${dir_path}/truvari/ \
      -f ${referenceFasta} \
      --includebed ${benchmark_bed} \
      --passonly \
      --pctsim=0 \
      -r 2000 --giabreport &> /dev/null
fi
rm -rf ${dir_path}/truvari_DEL/
if [ "$comp_against" = "DEL" ]; then
   truvari \
      -b ${benchmark_vcf_DEL} \
      -c ${dir_path}/${id}.vcf.gz \
      -o ${dir_path}/truvari_DEL/ \
      -f ${referenceFasta} \
      --includebed ${benchmark_bed} \
      --passonly \
      --pctsim=0 \
      -r 2000 --giabreport &> /dev/null
fi
rm -rf ${dir_path}/truvari_INS/
if [ "$comp_against" = "INS" ]; then
   truvari \
      -b ${benchmark_vcf_INS} \
      -c ${dir_path}/${id}.vcf.gz \
      -o ${dir_path}/truvari_INS/ \
      -f ${referenceFasta} \
      --includebed ${benchmark_bed} \
      --passonly \
      --pctsim=0 \
      -r 2000 --giabreport &> /dev/null
fi

echo "Done."

