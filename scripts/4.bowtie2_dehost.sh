# !/usr/bin/bash

############################################################################################################################################################
# Script ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/4.bowtie2_dehost.sh 

# Step 4 of the analysis.
# This script is for eliminating human reads contamination in the pre-processed metagenomic data with bowtie2 after filtering low-quality reads in script 2.
# Last modified: 23-6-13
############################################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/2.fastp
OUTPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
OUTPUT_BAM_DIR=${OUTPUT_DIR}/bam
OUTPUT_LOG_DIR=${OUTPUT_DIR}/log
BOWTIE2_REFERENCE=~/genome/human/bowtie2_reference/human

mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_BAM_DIR}
mkdir -p ${OUTPUT_LOG_DIR}

awk '{print $1"\t"$4"_"$5"_"$1}' ~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/download_scripts/sample_design.txt | grep "TD" | sed -e "s/\n//g" | sort | while read sample info
do
	cd ${INPUT_DIR}
        #echo ${sample} ${info}
        #rename "${sample}_" "${info}_" *
done

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi


# Mapping reads to human genome with bowtie2.
cd ${INPUT_DIR}
seq ${BEGIN} ${END}| while read individual
do
	echo ------------------------------------------------------------------------------
	echo ${individual}

	ls ${individual}_*_*.fastq.gz | sed -e 's/_fastp_1.fastq.gz//g' -e 's/_fastp_2.fastq.gz//g' -e 's/_unpaired.fastq.gz//g' | sort | uniq | while read sample
	do

		#input files
		PAIRED_1=${INPUT_DIR}/${sample}"_fastp_1.fastq.gz"
		PAIRED_2=${INPUT_DIR}/${sample}"_fastp_2.fastq.gz"
		UNPAIRED=${INPUT_DIR}/${sample}"_unpaired.fastq.gz"

		#output sam/bam files
		PE_SAM=${OUTPUT_BAM_DIR}/${sample}"_paired.sam"
		SE_SAM=${OUTPUT_BAM_DIR}/${sample}"_unpaired.sam"
		PE_BAM=${OUTPUT_BAM_DIR}/${sample}"_paired.bam"
        	SE_BAM=${OUTPUT_BAM_DIR}/${sample}"_unpaired.bam"
	
		#output unaligned reads
		PE_UNALIGNED=${OUTPUT_DIR}/${sample}"_paired_dehost.fastq.gz"
		SE_UNALIGNED=${OUTPUT_DIR}/${sample}"_unpaired_dehost.fastq.gz"
	
		#output log files
		PE_LOG=${OUTPUT_LOG_DIR}/${sample}"_paired.log"
		SE_LOG=${OUTPUT_LOG_DIR}/${sample}"_unpaired.log"
	
		echo ${sample}' paired'
		bowtie2 -p 16 --no-unal --un-conc-gz ${PE_UNALIGNED} -x ${BOWTIE2_REFERENCE} -1 ${PAIRED_1} -2 ${PAIRED_2} -S ${PE_SAM} 2> ${PE_LOG}
		#awk '($2==4) {print $1}' ${PE_SAM} > ${UNMAPPED_NAMES}
		#seqkit grep -f ${UNMAPPED_NAMES} ${FASTQ} > ${DEHOST_FASTQ}
		#gzip ${DEHOST_FASTQ}

		echo ${sample}' unpaired'
		bowtie2 -p 16 --no-unal --un-gz ${SE_UNALIGNED} -x ${BOWTIE2_REFERENCE} -U ${UNPAIRED} -S ${SE_SAM} 2> ${SE_LOG}

		#samtools view -b -S ${PE_SAM} > ${PE_BAM}
		#samtools view -b -S ${SE_SAM} > ${SE_BAM}
		rm ${PE_SAM}
		rm ${SE_SAM}
	done

done
