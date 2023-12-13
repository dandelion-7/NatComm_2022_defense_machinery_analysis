# !/usr/bin/bash

# Step 3 for illumina illumina sequencing analysis.
# This script is for quality control of illumina sequencing data after filtering with fastp.
# Last modified: 23-5-18.

# Set parameters.
BASE=~/crisprome/illumina_nanopore_comparing
INPUT_NANOPORE_DIR=~/public_raw_data/innerMongolian_2023_natMircobiol_Microbial_Genomes/raw/data/nanopore
INPUT_ILLUMINA_DIR=${BASE}/intermediates/illumina_fastp
OUTPUT_DIR=${BASE}/intermediates/illumina_fastqc_2
ADAPTOR_DIR=~/software/adaptor_sequences/adaptors.txt
mkdir -p ${OUTPUT_DIR}

ls ${INPUT_ILLUMINA_DIR} | sed '/unpaired/d' | sed -e 's/_fastp_1.fastq.gz//g' -e 's/_fastp_2.fastq.gz//g' | uniq | while read sample
do
	INPUT_FASTQ_1=${INPUT_ILLUMINA_DIR}/${sample}'_fastp_1.fastq.gz'
	INPUT_FASTQ_2=${INPUT_ILLUMINA_DIR}/${sample}'_fastp_2.fastq.gz'
	echo ${INPUT_FASTQ_1}
	echo ${INPUT_FASTQ_2}
	fastqc ${INPUT_FASTQ_1} ${INPUT_FASTQ_2} -o ${OUTPUT_DIR} -t 16 -a ${ADAPTOR_DIR}
done

# Perform multiqc
OUTPUT_MULTIQC_DIR=${OUTPUT_DIR}/multiqc
mkdir -p ${OUTPUT_MULTIQC_DIR}
multiqc -o ${OUTPUT_MULTIQC_DIR} ${OUTPUT_DIR}
