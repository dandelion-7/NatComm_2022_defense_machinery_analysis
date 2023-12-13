# !/usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script:~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/22.virsorter2.sh 
# This script is for predicting and extracting the viral contigs from the total assembled contigs in script 10.
# Last modified: 23-9-15.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/10.megahit
OUTPUT_DIR=${BASE}/intermediates/22.virsorter2
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate vs2
seq ${BEGIN} ${END} | while read individual
do
	echo -----------------------------------------------------------------------------------------
	echo ${individual}
	input=${INPUT_DIR}/${individual}/${individual}'_noted.contigs.fa'
	output=${OUTPUT_DIR}/${individual}
	virsorter run -w ${output} -i ${input} --include-groups "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae" --label ${individual} --keep-original-seq -j 64 --min-length 1000
done
conda deactivate
