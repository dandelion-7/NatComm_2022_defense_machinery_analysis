# !/usr/bin/bash

############################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/1.fastqc_1.sh
# Step 1.
# This script is for quality control of illumina sequencing data before filtering.
# Last modified: 23-6-13.
###########################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/data/illumina/time_series
OUTPUT_DIR=${BASE}/intermediates/1.fastqc_1
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
ls ${INPUT_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' | sort | uniq | sed -n "${BEGIN}, ${END}p" | while read sample
do
	INPUT_FASTQ_1=${INPUT_DIR}/${sample}'_1.fastq.gz'
	INPUT_FASTQ_2=${INPUT_DIR}/${sample}'_2.fastq.gz'
	#echo ${INPUT_FASTQ_1}
	#echo ${INPUT_FASTQ_2}
	fastqc ${INPUT_FASTQ_1} ${INPUT_FASTQ_2} -o ${OUTPUT_DIR} -t 64 -a ${ADAPTOR_DIR}
done

# Perform multiqc
OUTPUT_MULTIQC_DIR=${OUTPUT_DIR}/multiqc
#mkdir -p ${OUTPUT_MULTIQC_DIR}
#multiqc -o ${OUTPUT_MULTIQC_DIR} ${OUTPUT_DIR}
