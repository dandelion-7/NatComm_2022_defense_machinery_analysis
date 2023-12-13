# !/usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/25-4.merge_prodigal_results.sh
# This script is for merging the results of prodigal (predicted proteins and gens) for annotating the contigs.
# Last modified: 23.10.12
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/25.defensefinder/prodigal

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

touch ${INPUT_DIR}/total_CDS.txt
seq ${BEGIN} ${END} | while read individual
do
	echo ------------------------------------------------------------------------------------------------------------------------------------------------------------
	echo ${individual}
	cd ${INPUT_DIR}/${individual}_splitted_contigs
	ls *_proteins.fasta | sed 's/_proteins.fasta//g' | while read sample
	do
		echo ${sample}
		grep ">" ${INPUT_DIR}/${individual}_splitted_contigs/${sample}_proteins.fasta | sed -e 's/ //g' -e 's/#/;/g' -e 's/>//g' >> ${INPUT_DIR}/total_CDS.txt
	done
done
