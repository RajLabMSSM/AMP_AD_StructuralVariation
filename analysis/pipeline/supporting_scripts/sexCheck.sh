#!/bin/bash

#################################
#             HOLMES            #
#  The liWGS SV discovery tool  #
#################################

# Copyright (c) 2016 Ryan L. Collins and the laboratory of Michael E. Talkowski
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# Code development credits and citation availble on GitHub

module load sambamba

#Sex Check Script

#Read input
ID=$1
bam=$2
DICT=$3
runDir=$4

#Set additional params
genome_size=2897310462 #n-masked length of grch37

#Get primary alignment count in library
total=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} )

#Get primary alignment count to X and Y
X=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} X )
Y=$( sambamba view -c -F 'not secondary_alignment and not duplicate' ${bam} Y )

#Get expected fractions
# genome_size=$( fgrep -w "@SQ" $DICT | cut -f3 | cut -d\: -f2 | awk '{ sum+=$1 } END { print sum }' ) #dynamically calculate genome size, does not take into account n-masking
Xfrac_ex=$( echo -e "scale=6; (( $( fgrep X ${DICT} | cut -f3 | cut -d\: -f2 ) / ${genome_size} ))" | bc )
Yfrac_ex=$( echo -e "scale=6; (( $( fgrep Y ${DICT} | cut -f3 | cut -d\: -f2 ) / ${genome_size} ))" | bc )

#Get observed fractions
Xfrac_obs=$( echo -e "scale=6; (( ${X} / ${total} ))" | bc )
Yfrac_obs=$( echo -e "scale=6; (( ${Y} / ${total} ))" | bc )

#Predict sex
Xcopies=$( echo -e "scale=1; (( 2 * ${Xfrac_obs} / ${Xfrac_ex} ))" | bc | awk '{printf "%.0f\n",$1}' )
Ycopies=$( echo -e "scale=1; (( 2 * ${Yfrac_obs} / ${Yfrac_ex} ))" | bc | awk '{printf "%.0f\n",$1}' )
if [ ${Xcopies} -eq 2 ] && [ ${Ycopies} -eq 0 ]; then
  predSex="F"
elif [ ${Xcopies} -eq 1 ] && [ ${Ycopies} -eq 1 ]; then
  predSex="M"
else
  predSex="O"
fi

#Write Results to sexcheck file
echo -e "${ID} SEX CHECK RESULTS\n\nLibrary primary alignment count: ${total}\nX primary alignment count: ${X}\nY primary alignment count: ${Y}\
\n\nX library fraction observed: ${Xfrac_obs} (${Xfrac_ex} expected for diploid)\nY library fraction observed: ${Yfrac_obs} (${Yfrac_ex} expected for diploid)\
\n\nPredicted X copy state: ${Xcopies}\nPredicted Y copy state: ${Ycopies}\n\nPREDICTED SEX: ${predSex}" > ${runDir}/${ID}.sexCheck


