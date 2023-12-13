# !/usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/25-2.defensefinder.sh
# This script is for detecting anti-phage related genes from CDSs of assembled contigs of each individual with DefenseFinder following script 25-1.
# Last modified: 23.10.8
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/25.defensefinder/prodigal
OUTPUT_DIR=${BASE}/intermediates/25.defensefinder/defensefinder
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

head -n 1 ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/1/1_long_contigs.part_001/1_long_contigs.part_001_proteins_defense_finder_hmmer.tsv > ${OUTPUT_DIR}/merged_defensefinder_hmmer.tsv
head -n 1 ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/1/1_long_contigs.part_001/1_long_contigs.part_001_proteins_defense_finder_genes.tsv > ${OUTPUT_DIR}/merged_defensefinder_genes.tsv
head -n 1 ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/1/1_long_contigs.part_001/1_long_contigs.part_001_proteins_defense_finder_systems.tsv > ${OUTPUT_DIR}/merged_defensefinder_systems.tsv

source activate defensefinder
seq ${BEGIN} ${END} | while read individual
do
	echo ------------------------------------------------------------------------------------------------------------------------------------------------------------
	echo ${individual}
	cd ${INPUT_DIR}/${individual}_splitted_contigs
	ls *_proteins.fasta | sed 's/_proteins.fasta//g' | while read sample
	do
		echo ${sample}
		sed '1d' ${OUTPUT_DIR}/${individual}/${sample}/*_hmmer.tsv >> ${OUTPUT_DIR}/merged_defensefinder_hmmer.tsv
		sed '1d' ${OUTPUT_DIR}/${individual}/${sample}/*_genes.tsv >> ${OUTPUT_DIR}/merged_defensefinder_genes.tsv
		sed '1d' ${OUTPUT_DIR}/${individual}/${sample}/*_systems.tsv >> ${OUTPUT_DIR}/merged_defensefinder_systems.tsv
	done
done
conda deactivate
