# !/usr/bin/bash

# Step 1 for illumina illumina sequencing analysis.
# This script is for quality control of illumina sequencing data before filtering.
# Last modified: 23-5-17.

# Set parameters.
BASE=~/crisprome/illumina_nanopore_comparing
INPUT_NANOPORE_DIR=~/public_raw_data/innerMongolian_2023_natMircobiol_Microbial_Genomes/raw/data/nanopore
INPUT_ILLUMINA_DIR=~/public_raw_data/innerMongolian_2023_natMircobiol_Microbial_Genomes/raw/data/illumina
OUTPUT_DIR=${BASE}/intermediates/nanopore_fastqc_1
ADAPTOR_DIR=~/software/adaptor_sequences/adaptors.txt
mkdir -p ${OUTPUT_DIR}

ls ${INPUT_NANOPORE_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' | uniq | while read sample
do
	INPUT_FASTQ_1=${INPUT_NANOPORE_DIR}/${sample}'_1.fastq.gz'
	INPUT_FASTQ_2=${INPUT_NANOPORE_DIR}/${sample}'_2.fastq.gz'
	echo ${INPUT_FASTQ_1}
	echo ${INPUT_FASTQ_2}
	fastqc ${INPUT_FASTQ_1} ${INPUT_FASTQ_2} -o ${OUTPUT_DIR} -t 16
done

# Perform multiqc
OUTPUT_MULTIQC_DIR=${OUTPUT_DIR}/multiqc
mkdir -p ${OUTPUT_MULTIQC_DIR}
multiqc -o ${OUTPUT_MULTIQC_DIR} ${OUTPUT_DIR}
