# !usr/bin/bash
#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/14.checkm.sh
# This script is for assigning taxonomy of bins and evaluate the completeness/contamination of bins, following script 12.
# Last modified: 23-7-19.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/12.metabat2/bin
OUTPUT_DIR=${BASE}/intermediates/14.checkm
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate checkm
seq ${BEGIN} ${END} | while read individual
do 
	echo ${individual}
	checkm lineage_wf -t 32 --pplacer_threads 16 -x fa ${INPUT_DIR}/${individual} ${OUTPUT_DIR}
done
conda deactivate
