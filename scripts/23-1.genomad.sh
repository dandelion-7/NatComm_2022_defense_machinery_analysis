# !/usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script:~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/23-1.genomad.sh 
# This script is for predicting and extracting the viral contigs from the total assembled contigs following script 10, with geNomad.
# Last modified: 23-9-28.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/10.megahit
OUTPUT_DIR=${BASE}/intermediates/23.genomad/total_pipeline
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate genomad
seq ${BEGIN} ${END} | while read individual
do
	echo -----------------------------------------------------------------------------------------
	echo ${individual}
	input=${INPUT_DIR}/${individual}/${individual}'_noted.contigs.fa'
	output=${OUTPUT_DIR}/${individual}
	#genomad end-to-end --disable-nn-classification --cleanup --splits 8 ${input} ${output} ~/software/genomad/genomad_db
	genomad end-to-end --cleanup --splits 8 --threads 2 ${input} ${output} ~/software/genomad/genomad_db
done
conda deactivate
