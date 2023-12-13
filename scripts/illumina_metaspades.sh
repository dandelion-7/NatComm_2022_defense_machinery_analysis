#!/usr/bin/bash

# Step 4 of illumina sequencing data analysis.
# This script is for assembling the illumina shotgun metagenomic data after filtering with fastp with metaSPAdes.
# Last modified: 23-5-18.

# Set parameters.
BASE=~/crisprome/illumina_nanopore_comparing
INPUT_DIR=${BASE}/intermediates/illumina_fastp
OUTPUT_DIR=${BASE}/intermediates/illumina_metaspades
METASPADES=~/software/spades/SPAdes-3.15.5-Linux/bin/metaspades.py
mkdir -p ${OUTPUT_DIR}

# Perform assembly with SPAdes.
ls ${INPUT_DIR}| sed '/unpaired/d' | sed -e 's/_fastp_1.fastq.gz//g' -e 's/_fastp_2.fastq.gz//g' | uniq | while read sample
do	
	INPUT_PAIRED_1=${INPUT_DIR}/${sample}'_fastp_1.fastq.gz'
	INPUT_PAIRED_2=${INPUT_DIR}/${sample}'_fastp_2.fastq.gz'
	INPUT_UNPAIRED=${INPUT_DIR}/${sample}'_unpaired.fastq.gz'
	OUTPUT_SAMPLE_DIR=${OUTPUT_DIR}/${sample}
	mkdir -p ${OUTPUT_SAMPLE_DIR}
	echo ${INPUT_PAIRED_1}
	${METASPADES} -1 ${INPUT_PAIRED_1} -2 ${INPUT_PAIRED_2} -s ${INPUT_UNPAIRED} -t 48 --memory 150 -o ${OUTPUT_SAMPLE_DIR}
done

# RAM can't fulfill the requirements of metaspades, so error 12 is returned (a spades.log file is reserved in the OUTPUT_DIR). MEGAHIT will be used for assembling instead.
