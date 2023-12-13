# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/19.bowtie2_reads_contig.sh

# This script is for mapping the thoroguh reads corresponding to CRISPR from metagenomic reads onto the assembled contigsfor further assigning taxonomy to the CRISPR arrays.
# Last modified: 23-8-31.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_CONTIGS_DIR=${BASE}/intermediates/10.megahit
INPUT_READS_DIR=${BASE}/intermediates/5.crass
OUTPUT_DIR=${BASE}/intermediates/19.bowtie2_reads_contig
OUTPUT_NOTED_READS_DIR=${OUTPUT_DIR}/noted_reads
OUTPUT_REF_DIR=${OUTPUT_DIR}/references
OUTPUT_BAM_DIR=${OUTPUT_DIR}/bam
OUTPUT_LOG_DIR=${OUTPUT_DIR}/log
OUTPUT_UNALIGNED_DIR=${OUTPUT_DIR}/unaligned
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_NOTED_READS_DIR}
mkdir -p ${OUTPUT_REF_DIR}
mkdir -p ${OUTPUT_BAM_DIR}
mkdir -p ${OUTPUT_LOG_DIR}
mkdir -p ${OUTPUT_UNALIGNED_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

# Add the individual/timepoint/sample/CRISPR group information onto the name of each read corresponding to CRISPR arrays, which makes it easy to retrieve the belonging of each read after mapping.
awk '{print $1"\t"$4"_"$5"_"$1}' ~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/download_scripts/sample_design.txt | grep "TD" | sed -e "s/\n//g" | sort | while read sample info
do
	#echo ----------------------------------------------------
	#echo ${info}
	ls ${INPUT_READS_DIR}/${sample} | grep ".fa"  | while read i
	do
		group=${i%_*}
		cat ${INPUT_READS_DIR}/${sample}/${i} | sed -e "s/>/>${info}_${group}_/g" -e "s/ /_/g" > ${OUTPUT_NOTED_READS_DIR}/${info}_${i}
	done
done

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
	TOTAL=${OUTPUT_NOTED_READS_DIR}/${i}"_total_reads.fa"
	echo ----------------------------------------------------------
	echo ${i}
	#ls ${OUTPUT_NOTED_READS_DIR}/${i}_*_TD*.fa
	cat ${OUTPUT_NOTED_READS_DIR}/${i}_*_TD*.fa > ${TOTAL}
	#echo ${TOTAL}

	# 1st round, align CRISPR-corresponding reads to assembled contigs in end-to-end mode.
	SAM=${OUTPUT_BAM_DIR}/${i}"_total_reads.sam"
	BAM=${OUTPUT_BAM_DIR}/${i}"_total_reads.bam"
	SORTED_BAM=${OUTPUT_BAM_DIR}/${i}"_total_reads.sorted.bam"
	LOG=${OUTPUT_LOG_DIR}/${i}"_total_reads.log"
	UNALIGNED=${OUTPUT_UNALIGNED_DIR}/${i}"_unaligned_endTOend.fa"
	BED=${OUTPUT_BAM_DIR}/${i}"_total_reads.bed"
	bowtie2 -f -p 32 -x ${REFERENCE} -U ${TOTAL} --un ${UNALIGNED} -S ${SAM} 2> ${LOG}
	samtools view -b -S ${SAM} > ${BAM}
	#rm ${SAM}
	samtools sort --threads 32 -l 9 -o ${SORTED_BAM} ${BAM}
	bedtools bamtobed -i ${SORTED_BAM} > ${BED}
	rm ${TOTAL}

	# 2nd round, align the unaliged reads in round 1 with the local mode.
	SAM_2=${OUTPUT_BAM_DIR}/${i}"_unaligned_endTOend.sam"
        BAM_2=${OUTPUT_BAM_DIR}/${i}"_unaligned_endTOend.bam"
        SORTED_BAM_2=${OUTPUT_BAM_DIR}/${i}"_unaligned_endTOend.sorted.bam"
        LOG_2=${OUTPUT_LOG_DIR}/${i}"_unaligned_endTOend.log"
	UNALIGNED_2=${OUTPUT_UNALIGNED_DIR}/${i}"_unaligned_local.fa"
	BED_2=${OUTPUT_BAM_DIR}/${i}"_unaligned_endTOend.bed"
        bowtie2 --local -f -p 32 -x ${REFERENCE} -U ${UNALIGNED} --un ${UNALIGNED_2} -S ${SAM_2} 2> ${LOG_2}
        samtools view -b -S ${SAM_2} > ${BAM_2}
        #rm ${SAM_2}
        samtools sort --threads 32 -l 9 -o ${SORTED_BAM_2} ${BAM_2}
	bedtools bamtobed -i ${SORTED_BAM_2} > ${BED_2}
done
