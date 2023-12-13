# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/17-1.kraken_taxonomy.sh

# This script is for assigning taxonomy to the identified reads corresponding to CRISPR arrays with Kraken, following script 5.
# But results of this script are not used, because due to the presence of spacers, the taxonomic assignment with Kraken may not be reliable.
# Last modified: 23.8.17.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/5.crass
OUTPUT_DIR=${BASE}/intermediates/17.kraken_taxonomy/confidence_0.6
mkdir -p ${OUTPUT_DIR}	
KRAKEN=~/software/kraken2/kraken2_installation/kraken2
KRAKEN_DB=~/software/kraken2/kraken2_database

touch ${OUTPUT_DIR}/"total_unclassified.fasta"
touch ${OUTPUT_DIR}/"total_classified.fasta"
touch ${OUTPUT_DIR}/"total_output.txt"
touch ${OUTPUT_DIR}/"total_report.txt"

# Rename the outputs with the individual and time point numbers.
awk '{print $1"\t"$4"_"$5"_"$1}' ~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/download_scripts/sample_design.txt | grep "TD" | sed -e "s/\n//g" | sort | while read sample info
do
	OUTPUT_SAMPLE_DIR=${OUTPUT_DIR}/${info}
	mkdir -p ${OUTPUT_SAMPLE_DIR}
	INPUT_ARRAYS=${OUTPUT_SAMPLE_DIR}/${info}"_merged.fasta"
	touch ${INPUT_ARRAYS}

	echo '-----------------------------------------------------------------------------'

	#merge reads corresponding to CRISPR arrays from each sample into one fasta file.
	ls ${INPUT_DIR}/${sample} | grep "Group" | while read array
	do
		GROUP=${info}_${array%.*}
		cat ${INPUT_DIR}/${sample}/${array} | sed "s/>/>${GROUP}_/g" >> ${INPUT_ARRAYS}
	done


	UNCLASSIFIED_OUT=${OUTPUT_SAMPLE_DIR}/${info}"_unclassified.fasta"
        CLASSIFIED_OUT=${OUTPUT_SAMPLE_DIR}/${info}"_classified.fasta"
        OUTPUT=${OUTPUT_SAMPLE_DIR}/${info}"_output.txt"
	REPORT=${OUTPUT_SAMPLE_DIR}/${info}"_report.txt"
	echo ${sample}
	echo ${OUTPUT}
	${KRAKEN} --db ${KRAKEN_DB} --threads 64 --unclassified-out ${UNCLASSIFIED_OUT} --classified-out ${CLASSIFIED_OUT} --output ${OUTPUT} --report ${REPORT} --use-names ${INPUT_ARRAYS} --confidence 0.6

	cat ${UNCLASSIFIED_OUT} >> ${OUTPUT_DIR}/"total_unclassified.fasta"
	cat ${CLASSIFIED_OUT} >> ${OUTPUT_DIR}/"total_classified.fasta"
	cat ${OUTPUT} >> ${OUTPUT_DIR}/"total_output.txt"
	cat ${REPORT} >> ${OUTPUT_DIR}/"total_report.txt"
done
