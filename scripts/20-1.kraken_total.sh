# !/usr/bin/bash

#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/20-1.kraken_total.sh

# This script is for profiling the taxonomical composition of the whole metagenomic data, to help confirm the results of MMseqs2+bowtie2 taxonomy annotation.
# Last modified: 23.9.3.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
OUTPUT_DIR=${BASE}/intermediates/20.taxonomy_profiling/kraken
mkdir -p ${OUTPUT_DIR}
KRAKEN=~/software/kraken2/kraken2_installation/kraken2
KRAKEN_DB=~/software/kraken2/kraken2_database
echo -n > ${OUTPUT_DIR}/total_report.txt

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

# Rename the outputs with the individual and time point numbers.
seq ${BEGIN} ${END} | while read individual
do
	echo ------------------------------------------------------------------------
	echo ${individual}
	cd ${INPUT_DIR}
	ls ${individual}_*_*.fastq.gz | sed -e 's/_unpaired_dehost.fastq.gz//g' -e "s/_paired_dehost.1.fastq.gz//g" -e "s/_paired_dehost.2.fastq.gz//g" | uniq | while read sample
	do
		#echo ${sample}
		INPUT_1=${INPUT_DIR}/${sample}"_paired_dehost.1.fastq.gz"
		INPUT_2=${INPUT_DIR}/${sample}"_paired_dehost.2.fastq.gz"
		INPUT_UNPAIRED=${INPUT_DIR}/${sample}"_unpaired_dehost.fastq.gz"
		OUTPUT=${OUTPUT_DIR}/${sample}"_output.txt"
		REPORT=${OUTPUT_DIR}/${sample}"_report.txt" 
		${KRAKEN} --db ${KRAKEN_DB} --threads 32 --output ${OUTPUT} --report ${REPORT} --gzip-compressed ${INPUT_1} ${INPUT_2} ${INPUT_UNPAIRED} --confidence 0.6
		
		#sed 's/^[ ]*//g' ${REPORT}| sed 's/  //g' | sed 's/ /_/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | sed "s/$/&\t${sample}/g" | sed 's/ //g' >> ${OUTPUT_DIR}/total_report.txt
	done
done

cd ${OUTPUT_DIR}
ls *TD*_report.txt | sed 's/_report.txt//g' | while read sample
do
	echo ${sample}
	REPORT=${OUTPUT_DIR}/${sample}"_report.txt"
	sed 's/^[ ]*//g' ${REPORT}| sed 's/  //g' | sed 's/ /_/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | sed "s/$/&\t${sample}/g" | sed 's/ //g' >> ${OUTPUT_DIR}/total_report.txt
done
