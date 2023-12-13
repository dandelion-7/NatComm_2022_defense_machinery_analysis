# !/usr/bin/bash

######################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/cdhit.sh

# This script is for clustering and dereplicating the recovered DR and SP sequences with crisprtools from script 7.
# Last modified: 23-7-5.
######################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/7.crisprtools_extract
OUTPUT_DIR=${BASE}/intermediates/8.cdhit
mkdir -p ${OUTPUT_DIR}

echo "id	clstr	clstr_size	length	clstr_rep	clstr_iden	clstr_cov" > ${OUTPUT_DIR}/"repeats_clustered_summary.txt"
echo "id	clstr	clstr_size	length	clstr_rep	clstr_iden	clstr_cov" > ${OUTPUT_DIR}/"spacers_clustered_summary.txt"

source activate cdhit
seq 10 | while read i
do
	echo ${i}
	MERGED_REPEATS=${OUTPUT_DIR}/${i}"_repeats.fasta"
	MERGED_SPACERS=${OUTPUT_DIR}/${i}"_spacers.fasta"
	
	OUTPUT_REPEATS=${OUTPUT_DIR}/${i}"_repeats_clustered.fasta"
	OUTPUT_SPACERS=${OUTPUT_DIR}/${i}"_spacers_clustered.fasta"

	REPEATS_SUMMARY=${OUTPUT_DIR}/"repeats_clustered_summary.txt"
	SPACERS_SUMMARY=${OUTPUT_DIR}/"spacers_clustered_summary.txt"


	cat ${INPUT_DIR}/${i}"_"*"repeats.fa" > ${MERGED_REPEATS}
	cat ${INPUT_DIR}/${i}"_"*"spacers.fa" > ${MERGED_SPACERS}
	
	cd-hit -i ${MERGED_REPEATS} -o ${OUTPUT_REPEATS} -aL 0.95 -d 200 -c 0.95 -T 16 -M 16000
	cd-hit -i ${MERGED_SPACERS} -o ${OUTPUT_SPACERS} -aL 0.95 -d 200 -c 0.95 -T 16 -M 16000

	clstr2txt.pl ${OUTPUT_REPEATS}".clstr" | sed '/clstr_size/d' >> ${REPEATS_SUMMARY}
	clstr2txt.pl ${OUTPUT_SPACERS}".clstr" | sed '/clstr_size/d' >> ${SPACERS_SUMMARY}

done
conda deactivate
