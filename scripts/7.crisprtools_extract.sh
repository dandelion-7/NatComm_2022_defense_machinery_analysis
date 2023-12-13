# !/usr/bin/bash

###########################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/7.crisprtools_extract.sh

# Step 7.
# This script is for extracting spacer/repeats/flanker sequences from the CRISPR loci recovered by Crass, following script 5.
# Last modified: 23-6-21.
###########################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/5.crass
OUTPUT_DIR=${BASE}/intermediates/7.crisprtools_extract

mkdir -p ${OUTPUT_DIR}

source activate crass
awk '{print $1"\t"$4"_"$5"_"$1}' ~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/download_scripts/sample_design.txt | grep "TD" | sed -e "s/\n//g" | sort | while read sample info
do
	#echo ${sample}
	echo ${info}

	INPUT_SAMPLE_DIR=${INPUT_DIR}/${sample}
	CRISPR_FILE=${INPUT_SAMPLE_DIR}/crass.crispr
	
	CRISPR_STAT=${OUTPUT_DIR}/${info}_crispr_stats.txt
	SPACERS=${OUTPUT_DIR}/${info}_spacers.fa
	REPEATS=${OUTPUT_DIR}/${info}_repeats.fa
	FLANKERS=${OUTPUT_DIR}/${info}_flankers.fa
	#echo ${REPEATS}

	cd ${INPUT_SAMPLE_DIR}
	crisprtools stat -H ${CRISPR_FILE} > ${CRISPR_STAT}
	crisprtools extract -s ${CRISPR_FILE} | sed "s/>/>${info}_/g" > ${SPACERS}
	crisprtools extract -d ${CRISPR_FILE} | sed "s/>/>${info}_/g" > ${REPEATS}
	crisprtools extract -f ${CRISPR_FILE} | sed "s/>/>${info}_/g" > ${FLANKERS}
done

conda deactivate
