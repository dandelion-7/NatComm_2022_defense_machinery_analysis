# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/12.metabat2.sh

# This script is for binning the assembled contigs from each individual with the bam files through bowtie2 mapping, following script 10 (assembling contigs) and script 9(mapping with bowtie2)
# Last modified: 23-7-6.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_CONTIGS_DIR=${BASE}/intermediates/10.megahit
INPUT_BAM_DIR=${BASE}/intermediates/9.bowtie2/bam
OUTPUT_DIR=${BASE}/intermediates/12.metabat2
OUTPUT_DIR_BAM_DEPTH=${OUTPUT_DIR}/bam_depth
OUTPUT_DIR_BIN=${OUTPUT_DIR}/bin
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR_BAM_DEPTH}
mkdir -p ${OUTPUT_DIR_BIN}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate metabat2
seq ${BEGIN} ${END} | while read i
do
	depth=${OUTPUT_DIR_BAM_DEPTH}/${i}"_depth.txt"
	contig=${INPUT_CONTIGS_DIR}/${i}/${i}".contigs.fa"
	bin=${OUTPUT_DIR_BIN}/${i}
	mkdir -p ${bin}

	bam=`ls ${INPUT_BAM_DIR}/${i}_*_TD*.sorted.bam | sed 's/\n/ /g'`
	echo ${bam}
	echo ----------------------------------------------
	
	jgi_summarize_bam_contig_depths --minContigLength 500 --outputDepth ${depth} ${bam}
	metabat2 -m 1500 -i ${contig} -o ${bin}/${i} -a ${depth}
done

conda deactivate
