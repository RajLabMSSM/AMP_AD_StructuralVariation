#!/bin/bash

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

## SV Pipeline - 03_MergeTools
## Usage: ./run_03_MergeTools.sh myConfigFile.config

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

## This script create a consensus call merging different tools for each sample using best strategy found for each SV type
## Additionally, SVs overlapping centromeres and other regions are also removed (Bradler et al. 2016)

## Run folder
mkdir -p ${outputFolderPath} || exit 1
mkdir -p ${outputFolderPath}/03_MergeTools || exit 1
runDir=${outputFolderPath}"/03_MergeTools"
mkdir -p ${runDir} || exit 1
mkdir -p ${runDir}/merged || exit 1
mkdir -p ${runDir}/logs || exit 1

callsDir=${outputFolderPath}/02_SVtools

# =============================================================================
#                    Loop over each sample and merge
# =============================================================================

while [ 1 ]
do
    read bamFile || break
    fileName=$(basename "$bamFile")
    id=`echo $fileName | cut -d '.' -f 1`
    #echo "|-- $bamFile ---- $id --|"

    if [[ $(wc -l <${runDir}/logs/${id}/svBenchmark_AllTypes/ALL_merged.truvari) -ge 2 ]]; then
       echo "Rerunning? Cleaning folders"
       rm -rf ${runDir}/logs/${id}
    fi
    mkdir -p ${runDir}/logs/${id}

    # write a script for submiting each sample
    cat /dev/null > ${runDir}/logs/${id}/submit.sh
    chmod u+x ${runDir}/logs/${id}/submit.sh

    echo "ml bcftools" >> ${runDir}/logs/${id}/submit.sh
    echo "ml samtools" >> ${runDir}/logs/${id}/submit.sh
    echo "ml vcftools" >> ${runDir}/logs/${id}/submit.sh

    # ------------------------------------------------------ #
    # Merging DEL only... Genotyping info for DEL will be defined later.
    # A good combination was found using Manta and at least 2 other tools (75% recall and 92% precision).
    # So, we first select calls from Manta alone - overlapping calls within the Manta report are first merged using SURVIVOR.
    echo "echo \"Merging: DEL ------\"" >> ${runDir}/logs/${id}/submit.sh
    echo "find $callsDir -type f -name \"${id}.vcf.filt.DEL\" | grep \"svFiltered\" | grep 'Manta' > ${runDir}/logs/${id}/DEL_filtered.vcf.list1" >> ${runDir}/logs/${id}/submit.sh
    echo "${bench_truvari_script} ${runDir}/logs/${id}/DEL_filtered.vcf.list1 ${runDir}/logs/${id}/svBenchmark_DEL_1 1 DEL" >> ${runDir}/logs/${id}/submit.sh
    echo "cp ${runDir}/logs/${id}/svBenchmark_DEL_1/merged.vcf ${runDir}/logs/${id}/DEL_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh

    # Next we complement calls with other tools (at least 2 callers support)
    echo "find $callsDir -type f -name \"${id}.vcf.filt.DEL\" | grep \"svFiltered\" | grep -v 'Manta\|FusorSV' > ${runDir}/logs/${id}/DEL_filtered.vcf.list2" >> ${runDir}/logs/${id}/submit.sh
    echo "${bench_truvari_script} ${runDir}/logs/${id}/DEL_filtered.vcf.list2 ${runDir}/logs/${id}/svBenchmark_DEL_2 2 DEL" >> ${runDir}/logs/${id}/submit.sh 
    echo "cp ${runDir}/logs/${id}/svBenchmark_DEL_2/merged.vcf ${runDir}/logs/${id}/DEL_merged.vcf2" >> ${runDir}/logs/${id}/submit.sh

    # Merging all
    echo "find ${runDir}/logs/${id} -type f -name \"DEL_merged.vcf[1|2]\" > ${runDir}/logs/${id}/DEL_filtered.vcf.list" >> ${runDir}/logs/${id}/submit.sh
    echo "${bench_truvari_script} ${runDir}/logs/${id}/DEL_filtered.vcf.list ${runDir}/logs/${id}/svBenchmark_DEL 1 DEL" >> ${runDir}/logs/${id}/submit.sh

    # Fix genotypes and return only one column for this sample. Genotypes here will be replaced later using SV2!
    echo "cp ${runDir}/logs/${id}/svBenchmark_DEL/merged.vcf ${runDir}/logs/${id}/DEL_merged.vcf3" >> ${runDir}/logs/${id}/submit.sh
    echo "vcf-sort -c ${runDir}/logs/${id}/DEL_merged.vcf3 > ${runDir}/logs/${id}/DEL_merged.sorted.vcf3" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/DEL_merged.sorted.vcf3 > ${runDir}/logs/${id}/DEL_merged.sorted.vcf3.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "tabix -p vcf ${runDir}/logs/${id}/DEL_merged.sorted.vcf3.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "${svpipeline_r_lib}/fix_SURVIVORgenotypes.R ${runDir}/logs/${id}/DEL_merged.sorted.vcf3.gz ${id} ${runDir}/logs/${id}/DEL_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "gunzip -f ${runDir}/logs/${id}/DEL_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh

    # Copy final results into merged folder and remove temporaries
    echo "cp ${runDir}/logs/${id}/DEL_merged.vcf ${runDir}/merged/${id}.DEL_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    echo "rm ${runDir}/logs/${id}/DEL_merged.vcf1 ${runDir}/logs/${id}/DEL_merged.vcf2 ${runDir}/logs/${id}/DEL_merged.vcf3 ${runDir}/logs/${id}/DEL_merged.sorted.vcf3 ${runDir}/logs/${id}/DEL_merged.sorted.vcf3.gz ${runDir}/logs/${id}/DEL_merged.sorted.vcf3.gz.tbi" >> ${runDir}/logs/${id}/submit.sh

    # ------------------------------------------------------ #
    # Merging INS only
    # A simple UNION using Manta and BreakSeq calls (22% recall and 95% precision)
    echo "echo \"Merging: INS ------\"" >> ${runDir}/logs/${id}/submit.sh
    echo "find $callsDir -type f -name \"${id}.vcf.filt.INS\" | grep \"svFiltered\" | grep 'Manta\|BreakSeq' > ${runDir}/logs/${id}/INS_filtered.vcf.list" >> ${runDir}/logs/${id}/submit.sh
    # Merge tools keeping ALT?

    # Merge tools and benchmark
    echo "${bench_truvari_script} ${runDir}/logs/${id}/INS_filtered.vcf.list ${runDir}/logs/${id}/svBenchmark_INS 1 INS" >> ${runDir}/logs/${id}/submit.sh

    # Fix genotypes and return only one column for this sample
    echo "cp ${runDir}/logs/${id}/svBenchmark_INS/merged.vcf ${runDir}/logs/${id}/INS_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    echo "vcf-sort -c ${runDir}/logs/${id}/INS_merged.vcf1 > ${runDir}/logs/${id}/INS_merged.sorted.vcf1" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/INS_merged.sorted.vcf1 > ${runDir}/logs/${id}/INS_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "tabix -p vcf ${runDir}/logs/${id}/INS_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "${svpipeline_r_lib}/fix_SURVIVORgenotypes.R ${runDir}/logs/${id}/INS_merged.sorted.vcf1.gz ${id} ${runDir}/logs/${id}/INS_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "gunzip -f ${runDir}/logs/${id}/INS_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh

    # Copy final results into merged folder and remove temporaries
    echo "cp ${runDir}/logs/${id}/INS_merged.vcf ${runDir}/merged/${id}.INS_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    echo "rm ${runDir}/logs/${id}/INS_merged.vcf1 ${runDir}/logs/${id}/INS_merged.sorted.vcf1 ${runDir}/logs/${id}/INS_merged.sorted.vcf1.gz ${runDir}/logs/${id}/INS_merged.sorted.vcf1.gz.tbi" >> ${runDir}/logs/${id}/submit.sh

    # ------------------------------------------------------ #
    # Merging DUP calls
    # Delly, LUMPY, Manta, BreakDancer and CNVnator report DUPs. GIAB do not provide benchmark for DUP.
    # So, as a conservative approach, we keep calls found by at least 2 tools.
    echo "echo \"Merging: DUP ------\"" >> ${runDir}/logs/${id}/submit.sh
    echo "find $callsDir -type f -name \"${id}.vcf.filt.DUP\" | grep \"svFiltered\" | grep -v 'FusorSV' > ${runDir}/logs/${id}/DUP_filtered.vcf.list" >> ${runDir}/logs/${id}/submit.sh

    breakpoint_dist=1000 # max distance between breakpoints
    min_num_calls=2 # Minimum number of supporting caller
    use_type=1 # Take the type into account (1==yes, else no)
    use_strand=1 # Take the strands of SVs into account (1==yes, else no)
    dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
    min_sv_size=50 # Minimum size of SVs to be taken into account.

    echo "SURVIVOR merge ${runDir}/logs/${id}/DUP_filtered.vcf.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${filename} ${runDir}/logs/${id}/DUP_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    #echo "sed -i 's/FORMAT=<ID=DR,Number=1,Type=Integer/FORMAT=<ID=DR,Number=1,Type=String/g' ${runDir}/logs/${id}/DUP_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    #echo "sed -i 's/ID=LN,Number=1,Type=Integer/ID=LN,Number=1,Type=String,/g' ${runDir}/logs/${id}/DUP_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    #echo "sed -i 's/;AVGLEN=/;SVLEN=/g' ${runDir}/logs/${id}/DUP_merged.vcf" >> ${runDir}/logs/${id}/submit.sh

    # Fix genotypes and return only one column for this sample. Genotypes for DUP will be replaced by SV2 results.
    echo "vcf-sort -c ${runDir}/logs/${id}/DUP_merged.vcf1 > ${runDir}/logs/${id}/DUP_merged.sorted.vcf1" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/DUP_merged.sorted.vcf1 > ${runDir}/logs/${id}/DUP_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "tabix -p vcf ${runDir}/logs/${id}/DUP_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "${svpipeline_r_lib}/fix_SURVIVORgenotypes.R ${runDir}/logs/${id}/DUP_merged.sorted.vcf1.gz ${id} ${runDir}/logs/${id}/DUP_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "gunzip -f ${runDir}/logs/${id}/DUP_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh

    # Copy final results into merged folder and remove temporaries
    echo "cp ${runDir}/logs/${id}/DUP_merged.vcf ${runDir}/merged/${id}.DUP_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    echo "rm ${runDir}/logs/${id}/DUP_merged.vcf1 ${runDir}/logs/${id}/DUP_merged.sorted.vcf1 ${runDir}/logs/${id}/DUP_merged.sorted.vcf1.gz ${runDir}/logs/${id}/DUP_merged.sorted.vcf1.gz.tbi" >> ${runDir}/logs/${id}/submit.sh

    # ------------------------------------------------------ #
    # Merging INV calls
    # As for DUPs, GIAB do not provide benchmark for INV as well.
    # Delly, LUMPY, Manta and BreakDancer report INV. 
    # Therefore, we take at least 2 callers support for INV.
    echo "echo \"Merging: INV ------\"" >> ${runDir}/logs/${id}/submit.sh
    echo "find $callsDir -type f -name \"${id}.vcf.filt.INV\" | grep \"svFiltered\" | grep -v 'FusorSV' > ${runDir}/logs/${id}/INV_filtered.vcf.list" >> ${runDir}/logs/${id}/submit.sh

    breakpoint_dist=1000 # max distance between breakpoints
    min_num_calls=2 # Minimum number of supporting caller
    use_type=1 # Take the type into account (1==yes, else no)
    use_strand=1 # Take the strands of SVs into account (1==yes, else no)
    dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
    min_sv_size=50 # Minimum size of SVs to be taken into account.

    echo "SURVIVOR merge ${runDir}/logs/${id}/INV_filtered.vcf.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${filename} ${runDir}/logs/${id}/INV_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    #echo "sed -i 's/FORMAT=<ID=DR,Number=1,Type=Integer/FORMAT=<ID=DR,Number=1,Type=String/g' ${runDir}/logs/${id}/INV_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    #echo "sed -i 's/ID=LN,Number=1,Type=Integer/ID=LN,Number=1,Type=String,/g' ${runDir}/logs/${id}/INV_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    #echo "sed -i 's/;AVGLEN=/;SVLEN=/g' ${runDir}/logs/${id}/INV_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh # No need after SURVIVOR 1.0.6
    
    # Fix genotypes and return only one column for this sample
    echo "vcf-sort -c ${runDir}/logs/${id}/INV_merged.vcf1 > ${runDir}/logs/${id}/INV_merged.sorted.vcf1" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/INV_merged.sorted.vcf1 > ${runDir}/logs/${id}/INV_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "tabix -p vcf ${runDir}/logs/${id}/INV_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "${svpipeline_r_lib}/fix_SURVIVORgenotypes.R ${runDir}/logs/${id}/INV_merged.sorted.vcf1.gz ${id} ${runDir}/logs/${id}/INV_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "gunzip -f ${runDir}/logs/${id}/INV_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    
    # Copy final results into merged folder and remove temporaries
    echo "cp ${runDir}/logs/${id}/INV_merged.vcf ${runDir}/merged/${id}.INV_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    echo "rm ${runDir}/logs/${id}/INV_merged.vcf1 ${runDir}/logs/${id}/INV_merged.sorted.vcf1 ${runDir}/logs/${id}/INV_merged.sorted.vcf1.gz ${runDir}/logs/${id}/INV_merged.sorted.vcf1.gz.tbi" >> ${runDir}/logs/${id}/submit.sh

    # ------------------------------------------------------ #
    # Merging TRA/BND calls
    # No GIAB benchmarking as well.
    # TRAs is reported BreakDancer only! BND is reported by Manta, Lumpy and Delly.
    # At least 2 tools suport for our final callset.
    echo "echo \"Merging: TRA/BND ------\"" >> ${runDir}/logs/${id}/submit.sh
    echo "find $callsDir -type f -name \"${id}.vcf.filt.BND\" | grep \"svFiltered\" > ${runDir}/logs/${id}/TRA_BND_filtered.vcf.list" >> ${runDir}/logs/${id}/submit.sh
    echo "find $callsDir -type f -name \"${id}.vcf.filt.TRA\" | grep \"svFiltered\" >> ${runDir}/logs/${id}/TRA_BND_filtered.vcf.list" >> ${runDir}/logs/${id}/submit.sh

    breakpoint_dist=1000 # max distance between breakpoints
    min_num_calls=2 # Minimum number of supporting caller
    use_type=0 # Take the type into account (1==yes, else no) ------- IMPORTANT: SVTYPE is not considered here. TRA=BND?
    use_strand=1 # Take the strands of SVs into account (1==yes, else no)
    dist_based=0 # Estimate distance based on the size of SV (1==yes, else no).
    min_sv_size=0 # Minimum size of SVs to be taken into account. -------- Minimum SV SIZE changed to 0

    echo "SURVIVOR merge ${runDir}/logs/${id}/TRA_BND_filtered.vcf.list ${breakpoint_dist} ${min_num_calls} ${use_type} ${use_strand} ${dist_based} ${min_sv_size} ${filename} ${runDir}/logs/${id}/TRA_BND_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    # echo "sed -i 's/FORMAT=<ID=DR,Number=1,Type=Integer/FORMAT=<ID=DR,Number=1,Type=String/g' ${runDir}/logs/${id}/TRA_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    # echo "sed -i 's/ID=LN,Number=1,Type=Integer/ID=LN,Number=1,Type=String,/g' ${runDir}/logs/${id}/TRA_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    # echo "sed -i 's/;AVGLEN=/;SVLEN=/g' ${runDir}/logs/${id}/TRA_merged.vcf1" >> ${runDir}/logs/${id}/submit.sh
    
    # Fix genotypes and return only one column for this sample
    echo "vcf-sort -c ${runDir}/logs/${id}/TRA_BND_merged.vcf1 > ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1 > ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "tabix -p vcf ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "${svpipeline_r_lib}/fix_SURVIVORgenotypes.R ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1.gz ${id} ${runDir}/logs/${id}/TRA_BND_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "gunzip -f ${runDir}/logs/${id}/TRA_BND_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh

    # Copy final results into merged folder and remove temporaries
    echo "cp ${runDir}/logs/${id}/TRA_BND_merged.vcf ${runDir}/merged/${id}.TRA_BND_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    echo "rm ${runDir}/logs/${id}/TRA_BND_merged.vcf1 ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1 ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1.gz ${runDir}/logs/${id}/TRA_BND_merged.sorted.vcf1.gz.tbi" >> ${runDir}/logs/${id}/submit.sh

    # ------------------------------------------------------ #
    # Merging all SVTYPEs
    echo "echo \"Merging: ALL ------\"" >> ${runDir}/logs/${id}/submit.sh
    
    echo "bgzip -f -c ${runDir}/logs/${id}/DEL_merged.vcf > ${runDir}/logs/${id}/DEL_merged.vcf.gz; tabix -p vcf ${runDir}/logs/${id}/DEL_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/INS_merged.vcf > ${runDir}/logs/${id}/INS_merged.vcf.gz; tabix -p vcf ${runDir}/logs/${id}/INS_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/DUP_merged.vcf > ${runDir}/logs/${id}/DUP_merged.vcf.gz; tabix -p vcf ${runDir}/logs/${id}/DUP_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/INV_merged.vcf > ${runDir}/logs/${id}/INV_merged.vcf.gz; tabix -p vcf ${runDir}/logs/${id}/INV_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "bgzip -f -c ${runDir}/logs/${id}/TRA_BND_merged.vcf > ${runDir}/logs/${id}/TRA_BND_merged.vcf.gz; tabix -p vcf ${runDir}/logs/${id}/TRA_BND_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    
    echo "bcftools concat -a -O z -o ${runDir}/logs/${id}/ALL_merged.vcf.gz ${runDir}/logs/${id}/DEL_merged.vcf.gz ${runDir}/logs/${id}/INS_merged.vcf.gz ${runDir}/logs/${id}/DUP_merged.vcf.gz ${runDir}/logs/${id}/INV_merged.vcf.gz ${runDir}/logs/${id}/TRA_BND_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "tabix -p vcf ${runDir}/logs/${id}/ALL_merged.vcf.gz" >> ${runDir}/logs/${id}/submit.sh
    echo "zcat ${runDir}/logs/${id}/ALL_merged.vcf.gz > ${runDir}/logs/${id}/ALL_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    
    # Benchmark results
    echo "${bench_truvari} ${runDir}/logs/${id}/ALL_merged.vcf ${runDir}/logs/${id}/svBenchmark_AllTypes" >> ${runDir}/logs/${id}/submit.sh

    # Copy final file 
    echo "cp ${runDir}/logs/${id}/ALL_merged.vcf ${runDir}/merged/${id}.ALL_merged.vcf" >> ${runDir}/logs/${id}/submit.sh
    
    # Submit job
    if [[ ! $(wc -l <${runDir}/logs/${id}/svBenchmark_AllTypes/ALL_merged.truvari) -ge 2 ]]; then
       # failed. Rerun
       echo $id
       #${runDir}/logs/${id}/submit.sh
       bsub -n 1 -R "rusage[mem=6000]" -W 12:00 -oo ${runDir}/logs/${id}/job.out -eo ${runDir}/logs/${id}/job.err -P acc_ad-omics -q express -J ${id}.merge ${runDir}/logs/${id}/submit.sh
    fi   

done < $bamFileList

