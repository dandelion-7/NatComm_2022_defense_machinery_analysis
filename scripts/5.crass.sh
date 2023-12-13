# !/usr/bin/bash

##########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/5.crass.sh 

# Step 5.
# This script is for identifying CRISPR sequences from the unassembled illumina sequencing data after dehosting with script 4.
# Last modified: 23-6-18.
##########################################################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
OUTPUT_CRISPR_DIR=${BASE}/intermediates/5.crass
mkdir -p ${OUTPUT_CRISPR_DIR}


# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=`ls ${INPUT_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' | sort | uniq | wc -l`
else
	BEGIN=${1}
	END=${2}
fi


# Use Crass to identify CRISPR reads from the unassembled data.
source activate crass
ls ${INPUT_DIR} | sed -e '/unpaired/d' -e '/bam/d' -e '/log/d' | sed -e 's/_paired_dehost.1.fastq.gz//g' -e 's/_paired_dehost.2.fastq.gz//g' | uniq | sort | sed -n "${BEGIN}, ${END}p" | while read sample
do
	echo ${sample}
	PAIRED_1=${INPUT_DIR}/${sample}"_paired_dehost.1.fastq.gz"
	PAIRED_2=${INPUT_DIR}/${sample}"_paired_dehost.2.fastq.gz"
	UNPAIRED=${INPUT_DIR}/${sample}"_unpaired_dehost.fastq.gz"
	OUTPUT_SAMPLE_DIR=${OUTPUT_CRISPR_DIR}/${sample}
	mkdir ${OUTPUT_SAMPLE_DIR}

	crass --windowLength 6 --covCutoff 2 --kmerCount 6 --maxSpacer 100 --outDir ${OUTPUT_SAMPLE_DIR} ${PAIRED_1} ${PAIRED_2} ${UNPAIRED} -c red-blue --logLevel 4 -G -L
done
conda deactivate
