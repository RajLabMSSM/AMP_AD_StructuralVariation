#!/bin/bash

echo "sv type: $1"
echo "vcf file: $2"
echo "output: $3"
svtype=$1
vcfFile=$2
outputDir=$3
id_mapping=$4

# FastQTL
export PATH="~/ad-omics/ricardo/MyApps/fastqtl/bin:$PATH"

ml purge
ml R
ml samtools
ml bcftools
ml bedtools
ml vcftools

# Create output folder. All results and processing files will written inside
mkdir -p ${outputDir}

# Copy BED file with the expression info. Sort by position, compress and index
sortBed -header -i ${exprBedFile} > ${outputDir}/eQTLResidualExpression.bed; bgzip -c ${outputDir}/eQTLResidualExpression.bed > ${outputDir}/eQTLResidualExpression.bed.gz; tabix -p bed ${outputDir}/eQTLResidualExpression.bed.gz

# Add samples not present in VCF (if any) to exlusion list
bcftools query -l ${vcfFile} > ${outputDir}/samplesInVcf.list
head -n 1 ${outputDir}/eQTLResidualExpression.bed | cut --complement -f1,2,3,4 | sed 's/\t/\n/g' > ${outputDir}/samplesInBed.list
comm -3 <(sort ${outputDir}/samplesInVcf.list) <(sort ${outputDir}/samplesInBed.list) | sed 's/\t//g' > ${outputDir}/exclude_samples.exc
awk '!seen[$0]++' ${outputDir}/exclude_samples.exc > ${outputDir}/exclude_samples.uniq.exc

cat ${id_mapping} | cut -f1 | grep -v -x -f ${outputDir}/exclude_samples.uniq.exc > ${outputDir}/samplesToMerge.list
bcftools view --force-samples -S ${outputDir}/samplesToMerge.list ${vcfFile} | bcftools +fill-tags | bcftools sort -Oz -o ${outputDir}/${svtype}_merged_Final.vcf.gz

tabix -p vcf ${outputDir}/${svtype}_merged_Final.vcf.gz

# for fastQTL compatibility
ml gsl/2.2.1

# Run nominal Pass
for i in $(seq 1 22); do
fastQTL --vcf ${outputDir}/${svtype}_merged_Final.vcf.gz --bed ${outputDir}/eQTLResidualExpression.bed.gz --out ${outputDir}/${svtype}_fastQTL.nominal.output.chunk_$i.txt.gz \
        --chunk $i 22 \
        --window 1e6 \
        --seed 123456789 \
        --normal \
        --exclude-samples ${outputDir}/exclude_samples.uniq.exc
done
zcat ${outputDir}/${svtype}_fastQTL.nominal.output.chunk_*.txt.gz | gzip -c > ${outputDir}/${svtype}_FinalResults_fastQTL_nominal.gz
rm -rf ${outputDir}/${svtype}_fastQTL.nominal.output.chunk_*.txt.gz

# Run permutation pass
for i in $(seq 1 22); do
fastQTL --vcf ${outputDir}/${svtype}_merged_Final.vcf.gz --bed ${outputDir}/eQTLResidualExpression.bed.gz --out ${outputDir}/${svtype}_fastQTL.output.chunk_$i.txt.gz \
        --permute 10000 \
        --chunk $i 22 \
        --window 1e6 \
        --seed 123456789 \
        --normal \
        --exclude-samples ${outputDir}/exclude_samples.uniq.exc
done
zcat ${outputDir}/${svtype}_fastQTL.output.chunk_*.txt.gz | gzip -c > ${outputDir}/${svtype}_FinalResults_fastQTL.gz
rm -rf ${outputDir}/${svtype}_fastQTL.output.chunk_*.txt.gz

zcat ${outputDir}/${svtype}_FinalResults_fastQTL.gz > ${outputDir}/${svtype}_FinalResults_fastQTL.txt
