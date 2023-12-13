# !/usr/bin/bash

############################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/3.fastqc_2.sh
# Step 3.
# This script is for quality control of illumina sequencing data after filtering with fastp in script 2.
# Last modified: 23-6-13.
###########################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/2.fastp
OUTPUT_DIR=${BASE}/intermediates/3.fastqc_2
ADAPTOR_DIR=~/software/adaptor_sequences/adaptors.txt
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=`ls ${INPUT_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' | sort | uniq | wc -l`
else
	BEGIN=${1}
	END=${2}
fi

# Run fastqc.
ls ${INPUT_DIR} | sed -e 's/_fastp_1.fastq.gz//g' -e 's/_fastp_2.fastq.gz//g' -e 's/_unpaired.fastq.gz//g' | sort | uniq | sed -n "${BEGIN}, ${END}p" | while read sample
do
	INPUT_FASTQ_1=${INPUT_DIR}/${sample}'_fastp_1.fastq.gz'
	INPUT_FASTQ_2=${INPUT_DIR}/${sample}'_fastp_2.fastq.gz'
	INPUT_UNPAIRED=${INPUT_DIR}/${sample}'_unpaired.fastq.gz'
	#echo ${INPUT_FASTQ_1}
	#echo ${INPUT_FASTQ_2}
	fastqc ${INPUT_FASTQ_1} ${INPUT_FASTQ_2} ${INPUT_UNPAIRED} -o ${OUTPUT_DIR} -t 64 -a ${ADAPTOR_DIR}
done

# Perform multiqc
OUTPUT_MULTIQC_DIR=${OUTPUT_DIR}/multiqc
mkdir -p ${OUTPUT_MULTIQC_DIR}
multiqc -o ${OUTPUT_MULTIQC_DIR} ${OUTPUT_DIR}
