# !/usr/bin/bash

#####################################################################################################################################
# script:~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/10.megahit.sh

# This script is for assembling dehosted sequencing data with megahit, following script 4, parallele with script 6, because SPAdes runs too slow.
# Last modified: 23-7-2.
####################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
MERGED_DIR=${INPUT_DIR}/individual_merged
OUTPUT_DIR=${BASE}/intermediates/10.megahit
mkdir -p ${OUTPUT_DIR}
mkdir -p ${MERGED_DIR}

# Rename the sequencing files with the individual and time point numbers.
awk '{print $1"_\t"$4"_"$5"_"$1"_"}' ~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/download_scripts/sample_design.txt | grep "TD" | sed -e "s/\n//g" | sort | while read sample info
do
	cd ${INPUT_DIR}
        #echo ${sample} ${info}
        #rename "${sample}" "${info}" *
done


# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

# Concate each individual's data at ten timepoints, and run megahit to assemble contigs from the merged fastq files.
source activate megahit
seq ${BEGIN} ${END} | while read i
do
	MERGED_PE=${MERGED_DIR}/${i}"_merged_PE.fastq.gz"
	MERGED_SE=${MERGED_DIR}/${i}"_merged_SE.fastq.gz"
        echo ${i}
	ls ${INPUT_DIR}/${i}"_"*"_paired"*
	ls ${INPUT_DIR}/${i}"_"*"_unpaired"*
        #cat ${INPUT_DIR}/${i}"_"*"_paired"* > ${MERGED_PE}
	#cat ${INPUT_DIR}/${i}"_"*"_unpaired"* > ${MERGED_SE}

	echo "megahit -m 0.45 -t 32 --out-dir ${OUTPUT_DIR}/${i} --out-prefix ${i} --12 ${MERGED_PE} -r ${MERGED_SE}" #for megahit the output dir must be non-present before running.
	#megahit -m 0.45 -t 32 --out-dir ${OUTPUT_DIR}/${i} --out-prefix ${i} --12 ${MERGED_PE} -r ${MERGED_SE}
done
conda deactivate
