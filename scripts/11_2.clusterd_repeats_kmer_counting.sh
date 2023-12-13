# !/usr/bin/bash
############################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/11_2.clusterd_repeats_kmer_counting.sh

# This script is for calculating the k-mer profiles of stable repeats with jellyfish
# Last modified: 23-6-13.
###########################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/11.clustered_CRISPR_stats
OUTPUT_DIR=${BASE}/intermediates/11.clustered_CRISPR_stats/11_2.clustered_repeats_kmer_counting
mkdir -p ${OUTPUT_DIR}

source activate jellyfish
seq 1 10 | while read individual
do
	seq 1 10 | while read timepoint
	do
		echo "Individual:${individual}; Timepoint:${timepoint}"
		TEMP_INPUT_FASTA=${INPUT_DIR}/${individual}"_"${timepoint}"_stable_repeats.temp.fasta"
		grep ">${individual}_${timepoint}_TD" ${INPUT_DIR}/stable_repeats.fasta -A 1 > ${TEMP_INPUT_FASTA}

		for k in {4..15}
		do
			echo "${k}-mer"
			OUTPUT_JF=${OUTPUT_DIR}/${k}'_mer_'${individual}'_'${timepoint}'.jf'
			OUTPUT_FASTA=${OUTPUT_DIR}/${k}'_mer_'${individual}'_'${timepoint}'.fasta'
			OUTPUT_SUMMARY=${OUTPUT_DIR}/"stable_repeats_"${k}'_mer_summary.txt'

			jellyfish count -m ${k}	-s 100M -t 32 -C ${TEMP_INPUT_FASTA} -o ${OUTPUT_JF}
			jellyfish dump ${OUTPUT_JF} > ${OUTPUT_FASTA}
			grep ">" ${OUTPUT_FASTA} | sed 's/>//g' > ${OUTPUT_DIR}/"count.temp"
			grep -v ">" ${OUTPUT_FASTA} > ${OUTPUT_DIR}/"sequence.temp"
			paste ${OUTPUT_DIR}/"sequence.temp" ${OUTPUT_DIR}/"count.temp" > ${OUTPUT_DIR}/"seq_count.temp"
			awk '{print $1"\t"$2"\t'${individual}'\t'${timepoint}'"}' ${OUTPUT_DIR}/"seq_count.temp" >> ${OUTPUT_SUMMARY}

			rm ${OUTPUT_JF}
			rm ${OUTPUT_FASTA}
			rm ${OUTPUT_DIR}/*.temp
		done

		rm ${TEMP_INPUT_FASTA}
	done
done
