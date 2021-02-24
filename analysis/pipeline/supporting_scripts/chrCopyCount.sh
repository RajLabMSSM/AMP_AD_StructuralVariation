#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

#Full-Chromosome Aneuploidy Check

#Read input
ID=$1
bam=$2
DICT=$3
runDir=$4
NMASK=$5

#Set additional params
genome_size=2897310462 #n-masked length of grch37

#Write header to outfile
echo -e "#HOLMES ANEUPLOIDY CHECK\n#${ID}\n\n#contig\tlibrary_count\tcontig_count\t\
observed_fraction\texpected_fraction\tpredicted_copies" > \
${runDir}/${ID}.aneuploidyCheck

#Get primary alignment count in library
total=$( sambamba view -c -F 'proper_pair and not secondary_alignment and not duplicate' ${bam} )

#Iterate over primary contigs
for contig in $( seq 1 22 ) X Y; do
  length=$( echo "$( fgrep -w "SN:${contig}" ${DICT} | cut -f3 | \
  cut -d\: -f2 )-$( awk -v contig=${contig} '{ if ($1==contig) print $3-$2 }' ${NMASK} | \
  awk '{ sum+=$1 }END{ print sum }' )" | bc ) #get non-N-masked length of contig
  count=$( sambamba view -c -F 'proper_pair and not secondary_alignment and not duplicate' ${bam} ${contig} )
  frac_ex=$( echo -e "scale=6; (( ${length} / ${genome_size} ))" | bc )
  frac_obs=$( echo -e "scale=6; (( ${count} / ${total} ))" | bc )
  copies=$( echo -e "scale=3; (( 2 * ${frac_obs} / ${frac_ex} ))" | bc | \
    awk '{printf "%.2f\n",$1}' )
  echo -e "${contig}\t${total}\t${count}\t${frac_obs}\t${frac_ex}\t${copies}"
done >> ${runDir}/${ID}.aneuploidyCheck
