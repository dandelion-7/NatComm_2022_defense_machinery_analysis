# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/18-1.CRISPR_artifical_assembly.sh

# This script is for artifically assembling the recovered repeats and spacers by Crass from metagenomic data into an array, which will be suitble input for CRISPRidentify to evaluate the confidence level of this CRISPR array.
# Last modified: 23-8-23.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/7.crisprtools_extract
OUTPUT_DIR=${BASE}/intermediates/18.CRISPR_assembly_evaluate
OUTPUT_ARRAY_DIR=${BASE}/intermediates/18.CRISPR_assembly_evaluate/arrays
OUTPUT_EVALUATION_DIR=${BASE}/intermediates/18.CRISPR_assembly_evaluate/CRISPRidentify
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_ARRAY_DIR}
mkdir -p ${OUTPUT_EVALUATION_DIR}

# Artificially assemble/concate the identified CRISPR repeats/spacers into an "CRISPR array" as the input of evaluation by CRISPRidentify.
ls ${INPUT_DIR} | grep "repeats"| sed 's/_repeats.fa//g' | while read sample
do
	#echo ${sample}
	REPEATS=${INPUT_DIR}/${sample}_repeats.fa
	SPACERS=${INPUT_DIR}/${sample}_spacers.fa
	OUTPUT_ARRAY=${OUTPUT_ARRAY_DIR}/${sample}"_assembled_arrays.fa"
	#touch ${OUTPUT_ARRAY}

	grep ">" ${REPEATS} | sed 's/DR/\t/g' | awk '{print $1}' | while read group
	do
		repeat=`grep "${group}DR" -A 1 ${REPEATS} | grep -v ">"`
		spacer=`grep "${group}SP" -A 1 ${SPACERS} | grep -v ">" | tr "\n" "_"`
		sp_number=`grep "${group}SP" ${SPACERS} | wc -l`
		echo ${group}_${repeat}_${sp_number} >> ${OUTPUT_ARRAY}
		echo -n ${repeat} >> ${OUTPUT_ARRAY}; echo ${spacer} | sed "s/_/${repeat}/g" >> ${OUTPUT_ARRAY}
	done

done

#########################################################################################################################################################################

# Evaluate all the artificially assembled CRISPR arrays in one folder OUTPUT_ARRAY_DIR with the --input_folder_multifasta input mode of CRISPRidentify.
source activate crispr-softwares
python ~/software/CRISPRIdentify/CRISPRIdentify/CRISPRidentify-1.2.1/CRISPRidentify.py --input_folder_multifasta ${OUTPUT_ARRAY_DIR} --model 8 --result_folder ${OUTPUT_EVALUATION_DIR}
conda deactivate

#########################################################################################################################################################################

# Evaluate each individual's artificially assembled arrays separately with the --file input mode of CRISPRidentify.

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
        BEGIN=1
        END=10
else
        BEGIN=${1}
        END=${2}
fi

source activate crispr-softwares
seq ${BEGIN} ${END} | while read i
do
        INDIVIDUAL_FASTA=${OUTPUT_ARRAY_DIR}/${i}"_total_arrays.fasta"
        INDIVIDUAL_OUT_DIR=${OUTPUT_EVALUATION_DIR}/${i}"_model_10"
        mkdir -p ${INDIVIDUAL_OUT_DIR}
        echo --------------------------------------------------------
        echo ${i}
        #ls ${OUTPUT_ARRAY_DIR}/${i}_*_*.fa
        cat ${OUTPUT_ARRAY_DIR}/${i}_*_*.fa > ${INDIVIDUAL_FASTA}
        python ~/software/CRISPRIdentify/CRISPRIdentify/CRISPRidentify-1.2.1/CRISPRidentify.py --file ${INDIVIDUAL_FASTA} --model 10 --result_folder ${INDIVIDUAL_OUT_DIR}
done
conda deactivate
