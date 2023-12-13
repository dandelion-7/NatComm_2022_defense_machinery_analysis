# !/usr/bin/bash
#---------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/26-3.merge_featureCounts_results.sh
# This script is for merging the featureCount output matrix for downstream summary, following script 26-2.
# Last modified: 23.10.15.
#---------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/26.antivirus_genes_analysis/featureCounts_new

cd ${INPUT_DIR}
cat /dev/null > ${INPUT_DIR}/merged_gene_count.txt
#touch ${INPUT_DIR}/merged_gene_count.txt
seq 1 10 | while read individual
do
	echo ${individual}
	matrix=${INPUT_DIR}/${individual}_gene_count.txt
	sed '1,2d' ${matrix} | awk -v OFS="\t" '{print $1,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,'"${individual}"'}' | sed -e 's/+/\t/g' -e 's/:/=/g' >> ${INPUT_DIR}/merged_gene_count.txt
done
