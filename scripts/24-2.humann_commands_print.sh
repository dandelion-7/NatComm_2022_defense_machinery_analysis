# !/usr/bin/bash
#---------------------------------------------------------------------------------------------------
# Script:~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/24-1.humann.sh
# This script is for detecting Cas proteins from metagenome data with humann according to Uniref90 database.
# Last modified: 23.10.5
#---------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
OUTPUT_DIR=${BASE}/intermediates/24.cas_profiling/humann
mkdir -p ${OUTPUT_DIR}
BOWTIE2DB=~/software/metaphlan/bowtie2db_for_humann

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate humann
seq ${BEGIN} ${END} | while read individual
do
	#echo ------------------------------------------------------------------------
	#echo ${individual}
	cd ${INPUT_DIR}
	ls ${individual}_*_*.fastq.gz | sed -e 's/_unpaired_dehost.fastq.gz//g' -e "s/_paired_dehost.1.fastq.gz//g" -e "s/_paired_dehost.2.fastq.gz//g" | uniq | while read sample
	do
		echo ----------------------------------------------------------------------
		echo ${sample}
		echo gunzip ${INPUT_DIR}/${sample}"_paired_dehost.1.fastq.gz" -k
                echo gunzip ${INPUT_DIR}/${sample}"_paired_dehost.2.fastq.gz" -k
                echo gunzip ${INPUT_DIR}/${sample}"_unpaired_dehost.fastq.gz" -k
		
		INPUT_1=${INPUT_DIR}/${sample}"_paired_dehost.1.fastq"
		INPUT_2=${INPUT_DIR}/${sample}"_paired_dehost.2.fastq"
		INPUT_UNPAIRED=${INPUT_DIR}/${sample}"_unpaired_dehost.fastq"
		INPUT_MERGED=${OUTPUT_DIR}/${sample}"_merged.fastq"
		echo cat ${INPUT_1} ${INPUT_2} ${INPUT_UNPAIRED} > ${INPUT_MERGED} 
		echo rm ${INPUT_1}; echo rm ${INPUT_2}; echo rm ${INPUT_UNPAIRED}

		echo cd ${OUTPUT_DIR}
		echo humann -r --input ${INPUT_MERGED} --output-basename ${sample} --output ${OUTPUT_DIR}/${sample} --threads 16 --memory-use maximum --metaphlan-options "--bowtie2db /home/zhanggaopu/software/metaphlan/bowtie2db_for_humann --index mpa_vJan21_CHOCOPhlAnSGB_202103 --nproc 16" --bowtie-options "--threads 16" --diamond-options "--threads 16" --nucleotide-database ~/software/humann/humann_database/chocophlan --protein-database ~/software/humann/humann_database/uniref90/uniref
		echo rm ${INPUT_MERGED}
	done
done

conda deactivate
