# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/16.bowtie2_repeats_contig.sh

# This script is for mapping the recovered repeats from metagenomic reads onto the assembled contigs for further assigning taxonomy to the CRISPR arrays.
# But results of this script are not used, because the thorough CRISPR reads are also aligned to contigs in script 19.
# Last modified: 23-8-12.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_CONTIGS_DIR=${BASE}/intermediates/10.megahit
INPUT_REPEATS_DIR=${BASE}/intermediates/7.crisprtools_extract
OUTPUT_DIR=${BASE}/intermediates/16.bowtie2_repeats_contig
OUTPUT_REF_DIR=${OUTPUT_DIR}/references
OUTPUT_BAM_DIR=${OUTPUT_DIR}/bam
OUTPUT_LOG_DIR=${OUTPUT_DIR}/log
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_REF_DIR}
mkdir -p ${OUTPUT_BAM_DIR}
mkdir -p ${OUTPUT_LOG_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

# Construct reference for bowtie2
seq ${BEGIN} ${END} | while read i
do
	INPUT_CONTIGS=${INPUT_CONTIGS_DIR}/${i}/${i}".contigs.fa"
	REFERENCE=${i}"_contigs"

	cd ${OUTPUT_REF_DIR}
	#bowtie2-build -f ${INPUT_CONTIGS} ${REFERENCE}
done

# Bowtie2 mapping
seq ${BEGIN} ${END} | while read i
do
	REFERENCE=${OUTPUT_REF_DIR}/${i}"_contigs"
	TOTAL=${INPUT_REPEATS_DIR}/${i}"_total_repeats.fa"

	#ls ${INPUT_REPEATS_DIR}/${i}_*_TD*_repeats.fa
	cat ${INPUT_REPEATS_DIR}/${i}_*_TD*_repeats.fa > ${TOTAL}
	echo ${TOTAL}

	SAM=${OUTPUT_BAM_DIR}/${i}"_total_repeats.sam"
	BAM=${OUTPUT_BAM_DIR}/${i}"_total_repeats.bam"
	SORTED_BAM=${OUTPUT_BAM_DIR}/${i}"_total_repeats.sorted.bam"
	LOG=${OUTPUT_LOG_DIR}/${i}"_total_repeats.log"
	bowtie2 -f -p 32 -x ${REFERENCE} -U ${TOTAL} -S ${SAM} 2> ${LOG}
	samtools view -b -S ${SAM} > ${BAM}
	rm ${SAM}
	samtools sort --threads 32 -l 9 -o ${SORTED_BAM} ${BAM}
	rm ${TOTAL}

done
