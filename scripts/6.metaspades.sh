# !/usr/bin/bash

#####################################################################################################################################
# script:~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/6.metaspades.sh

# Step 6.
# This script is for assembling dehosted virome sequencing data with metaSPAdes, following script 4. Because MegaHit is used for mixed assembly, so the results of SPAdes were deleted.
# Last modified: 23-6-21.
####################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/4.bowtie2_dehost
OUTPUT_DIR=${BASE}/intermediates/6.metaspades
mkdir -p ${OUTPUT_DIR}

SPADES=~/software/spades/SPAdes-3.15.5-Linux/bin/spades.py

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=`ls ${INPUT_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' | sort | uniq | wc -l`
else
	BEGIN=${1}
	END=${2}
fi


# Perform assembly with SPAdes on each group of each individual.
ls ${INPUT_DIR} | grep "_paired_dehost.1.fastq.gz" | sed "s/_paired_dehost.1.fastq.gz//g" | sort | uniq | sed -n "${BEGIN}, ${END}p" | while read sample
do
	echo ${sample}
	# spades only accepts input files with the fastq/fastq.gz extensions, no numbers should be interwined.
	PAIRED_1=${INPUT_DIR}/${sample}'_paired_dehost.1.fastq.gz'
	PAIRED_2=${INPUT_DIR}/${sample}'_paired_dehost.2.fastq.gz'
	UNPAIRED=${INPUT_DIR}/${sample}'_unpaired_dehost.fastq.gz'
	OUTPUT_SAMPLE_DIR=${OUTPUT_DIR}/${sample}
	#mkdir -p ${OUTPUT_SAMPLE_DIR}

	#python ${SPADES} -1 ${PAIRED_1} -2 ${PAIRED_2} -s ${UNPAIRED} -o ${OUTPUT_SAMPLE_DIR} --meta --threads 75 --memory 155
done

# Rename the outputs with the individual and time point numbers.
awk '{print $1"\t"$4"_"$5"_"$1}' ~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/download_scripts/sample_design.txt | grep "TD" | sed -e "s/\n//g" | sort | while read sample info
do
	echo ${sample} ${info}
	#mv ${OUTPUT_DIR}/${sample}/contigs.fasta ${OUTPUT_DIR}/${sample}/${info}"_contigs.fasta"
	#mv ${OUTPUT_DIR}/${sample} ${OUTPUT_DIR}/${info}
	sed -i "s/>/>${info}_/g" ${OUTPUT_DIR}/${info}/${info}"_contigs.fasta"
done
