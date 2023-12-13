# !/usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/25-1.prodigal.sh
# This script is for predicting all the CDS from assembled contigs following script 10, for defense-finder.
# Last modified: 23.10.8
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/10.megahit
OUTPUT_DIR=${BASE}/intermediates/25.defensefinder/prodigal
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate defensefinder
seq ${BEGIN} ${END} | while read individual
do
	echo -----------------------------------------------------------------------------------------------
	echo ${individual}
	cd ${OUTPUT_DIR}
	contigs=${INPUT_DIR}/${individual}/${individual}_noted.contigs.fa
	long_contigs=${OUTPUT_DIR}/${individual}_long_contigs.fa
	
	echo "Extracting >500bp contigs with seqkit seq."
	seqkit seq -g -m 500 ${contigs} > ${long_contigs}

	echo "Split all the long contigs with seqkit split2."
	seqkit split2 ${long_contigs} --by-part 500 --out-dir ${OUTPUT_DIR}/${individual}_splitted_contigs -f -w 0 -j 16 --quiet

	echo "Predict CDSs with prodigal."
	cd ${OUTPUT_DIR}/${individual}_splitted_contigs
	ls *.fa | sed 's/.fa//g' | while read sample
	do
		echo ${sample}
		prodigal -i ${sample}.fa -a ${sample}_proteins.fasta -d ${sample}_genes.fasta -o ${sample}_output.txt -s ${sample}_potential_genes.txt -p meta -q
	done

done
conda deactivate
