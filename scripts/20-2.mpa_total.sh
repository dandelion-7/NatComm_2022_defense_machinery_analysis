# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/20-2.mpa_total.sh

# This script is for profiling the taxonomical composition of the whole metagenomic data with Metaphlan4, to help confirm the results of MMseqs2+bowtie2 taxonomy annotation.
# Last modified: 23.9.3.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
OUTPUT_DIR=${BASE}/intermediates/20.taxonomy_profiling/metaphlan
OUTPUT_BOWTIE2_DIR=${OUTPUT_DIR}/bowtie2out
OUTPUT_MPA_DIR=${OUTPUT_DIR}/mpa
mkdir -p ${OUTPUT_MPA_DIR}
mkdir -p ${OUTPUT_BOWTIE2_DIR}
BOWTIE2_DB=~/software/metaphlan/bowtie2db

echo -n > ${OUTPUT_MPA_DIR}/total_output.txt

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate mpa
seq ${BEGIN} ${END} | while read individual
do
	echo ------------------------------------------------------------------------
	echo ${individual}
	cd ${INPUT_DIR}
	ls ${individual}_*_*.fastq.gz | sed -e 's/_unpaired_dehost.fastq.gz//g' -e "s/_paired_dehost.1.fastq.gz//g" -e "s/_paired_dehost.2.fastq.gz//g" | uniq | while read sample
	do
		INPUT_1=${INPUT_DIR}/${sample}"_paired_dehost.1.fastq.gz"
		INPUT_2=${INPUT_DIR}/${sample}"_paired_dehost.2.fastq.gz"
		INPUT_UNPAIRED=${INPUT_DIR}/${sample}"_unpaired_dehost.fastq.gz"
		OUTPUT=${OUTPUT_MPA_DIR}/${sample}"_output.txt"
		BOWTIE2OUT=${OUTPUT_BOWTIE2_DIR}/${sample}"_bowtie2out.txt"
		echo ${sample}
		echo ${OUTPUT}

		cd ${OUTPUT_MPA_DIR}
		#metaphlan ${INPUT_1},${INPUT_2} --input_type fastq --bowtie2out ${BOWTIE2OUT} -o ./${sample}_output.txt --bowtie2db ${BOWTIE2_DB} --nproc 32 
		# If multiple files are input for metaphlan, they need to be sparated with ",", otherwise the following file will be re-written as the output.
		cat ./${sample}_output.txt | sed '1,5d' | awk '{print $1"\t"$2"\t"$3}' | sed -e 's/|/,/g' -e 's/ //g' | sed "s/$/&\t${sample}/g" >> ${OUTPUT_MPA_DIR}/total_output.txt
	done
done
