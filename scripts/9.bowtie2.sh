# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/9.bowtie2.sh

# This script is for mapping the filtered reads back to the contigs assembled with megahit from script 10, so as to get the relative abundancies of contigs for further binning and taxonomic analysis of the whole microbiome and recovered CRISPR reads.
# Last modified: 23-7-1.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/10.megahit
OUTPUT_DIR=${BASE}/intermediates/9.bowtie2
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

seq ${BEGIN} ${END} | while read i
do
	INPUT_CONTIGS=${INPUT_DIR}/${i}/${i}".contigs.fa"
	LONG_CONTIGS=${OUTPUT_REF_DIR}/${i}"_long_contigs.fasta"
	SHORT_CONTIGS=${OUTPUT_REF_DIR}/${i}"_short_contigs.fasta"
	REFERENCE=${i}"_long_contigs"
	
	#source activate seqkit
	#seqkit seq -g -m 500 ${INPUT_CONTIGS} > ${LONG_CONTIGS}
	#seqkit seq -g -M 499 ${INPUT_CONTIGS} > ${SHORT_CONTIGS}
	#conda deactivate

	#cd ${OUTPUT_REF_DIR}
	#bowtie2-build -f ${LONG_CONTIGS} ${REFERENCE}
done

INPUT_FASTA_DIR=${BASE}/intermediates/4.bowtie2_dehost

seq ${BEGIN} ${END} | while read i
do
	REFERENCE=${OUTPUT_REF_DIR}/${i}"_long_contigs"

	ls ${INPUT_FASTA_DIR}/${i}_*_TD* | sed -e "s#${INPUT_FASTA_DIR}/##g" -e 's/_paired_dehost.1.fastq.gz//g' -e 's/_paired_dehost.2.fastq.gz//g' -e 's/_unpaired_dehost.fastq.gz//g' | uniq | sort | while read sample
do
	PE_1=${INPUT_FASTA_DIR}/${sample}"_paired_dehost.1.fastq.gz"
	PE_2=${INPUT_FASTA_DIR}/${sample}"_paired_dehost.2.fastq.gz"
	SE=${INPUT_FASTA_DIR}/${sample}"_unpaired_dehost.fastq.gz"
	SAM=${OUTPUT_BAM_DIR}/${sample}".sam"
	BAM=${OUTPUT_BAM_DIR}/${sample}".bam"
	SORTED_BAM=${OUTPUT_BAM_DIR}/${sample}".sorted.bam"
	LOG=${OUTPUT_LOG_DIR}/${sample}".log"
	STATS=${OUTPUT_BAM_DIR}/${sample}"_stats.txt"
	echo ${sample}
	#bowtie2 -p 32 -x ${REFERENCE} -1 ${PE_1} -2 ${PE_2} -S ${SAM} 2> ${LOG}
	#samtools view -b -S ${SAM} > ${BAM}
	#rm ${SAM}
	#samtools sort --threads 32 -l 9 -o ${SORTED_BAM} ${BAM}
	#rm ${BAM}
	#samtools view --threads 32 ${SORTED_BAM} | awk '{if ($5 > 10) print}' | awk '{print $3}' | uniq -c | sed "s/$/&\t${sample}/g" > ${STATS}
	sed -i -e 's/ //g' -e 's/k/\tk/g' ${STATS}
done
	cd ${OUTPUT_BAM_DIR}
	#cat ${i}_*_stats.txt | sed -e 's/ //g' -e 's/k/\tk/g' > ${i}_stats.txt

done

cd ${OUTPUT_BAM_DIR}
#cat *TD*.txt | sed -e 's/ //g' -e 's/k/\tk/g' > total_stats.txt
