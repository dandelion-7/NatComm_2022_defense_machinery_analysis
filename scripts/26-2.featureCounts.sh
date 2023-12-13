# !/usr/bin/bash
#---------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/26.2-featureCounts.sh
# This script is for counting alignment over annotated antivirus genes with featureCounts following script 9 and 26-1.
# Last modified: 23.10.14.
#---------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
BAM_INPUT=${BASE}/intermediates/9.bowtie2/bam
GFF_INPUT=${BASE}/intermediates/26.antivirus_genes_analysis/gff_new
OUTPUT_DIR=${BASE}/intermediates/26.antivirus_genes_analysis/featureCounts_new; mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

seq ${BEGIN} ${END} | while read individual
do
	echo -----------------------------------------------------------------------------------------------------
	echo ${individual}
	#gff=${GFF_INPUT}/26-1.${individual}.gff # gff files from the old script 26-1.
	gff=${GFF_INPUT}/${individual}.gff # new gff files were generated with a new script 26-1.
	
	cd ${BAM_INPUT}
	samples=`ls ${individual}_*_*.sorted.bam`
	tmp=${OUTPUT_DIR}/${individual}_temp; mkdir -p ${tmp}
	echo ${samples}
	featureCounts -a ${gff} --tmpDir ${tmp} -t CDS -g gene_name -O -M -T 32 -p ${samples} -o ${OUTPUT_DIR}/${individual}_gene_count.txt
	rm -r ${tmp}
done
