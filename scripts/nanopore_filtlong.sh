# !/usr/bin/bash

# Step 1 for nanopore sequencing analysis.
# This script is for filtering low-quality and short reads in the nanopore sequencing file with filtlong.
# Last modified: 23-05-17.

# Set parameters.
BASE=~/crisprome/illumina_nanopore_comparing
INPUT_NANOPORE_DIR=~/public_raw_data/innerMongolian_2023_natMircobiol_Microbial_Genomes/raw/data/nanopore
INPUT_ILLUMINA_DIR=${BASE}/intermediates/illumina_fastp
OUTPUT_DIR=${BASE}/intermediates/nanopore_filtlong
mkdir -p ${OUTPUT_DIR}

source activate filtlong
ls ${INPUT_NANOPORE_DIR}|sed 's/_1.fastq.gz//g' | while read sample
do
	INPUT_ILLUMINA_FASTQ_1=${INPUT_ILLUMINA_DIR}/'Sample_'${sample}'_fastp_1.fastq.gz'
	INPUT_ILLUMINA_FASTQ_2=${INPUT_ILLUMINA_DIR}/'Sample_'${sample}'_fastp_2.fastq.gz'
	INPUT_NANOPORE_FASTQ=${INPUT_NANOPORE_DIR}/${sample}'_1.fastq.gz'
	OUTPUT_NANOPORE_FASTQ=${OUTPUT_DIR}/${sample}'_filtlong_1.fastq.gz'
	#echo ${INPUT_NANOPORE_FASTQ}
	#echo ${INPUT_ILLUMINA_FASTQ_1}
	#filtlong -1 ${INPUT_ILLUMINA_FASTQ_1} -2 ${INPUT_ILLUMINA_FASTQ_2} --min_length 1000 --mean_q_weight 6 --trim ${INPUT_NANOPORE_FASTQ} | gzip > ${OUTPUT_NANOPORE_FASTQ}
	filtlong --min_length 1000 --mean_q_weight 6 ${INPUT_NANOPORE_FASTQ} | gzip > ${OUTPUT_NANOPORE_FASTQ}
done
conda deactivate
