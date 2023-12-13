# !/usr/bin/bash 

#############################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/2.fastp.sh
# Step 2.
# This script is for filtering low-quality reads and removing adaptors in illumina sequencing data with fastp.
# Last modified: 23-6-13.

# The output files were deleted after finishing script 4.bowtie2_dehost.sh for saving the space. Re-run this script can get the output files again.
#############################################################################################################################################################

# Set parameters.
BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/data/illumina/time_series
OUTPUT_DIR=${BASE}/intermediates/2.fastp
ADAPTOR_DIR=~/software/adaptor_sequences/adaptors.txt
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=`ls ${INPUT_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' | sort | uniq | wc -l`
else
	BEGIN=${1}
	END=${2}
fi

echo ${BEGIN}
echo ${END}

source activate fastp
ls ${INPUT_DIR} | sed -e 's/_1.fastq.gz//g' -e 's/_2.fastq.gz//g' -e 's/META20IMWJ27_//g' | sort | uniq | sed -n "${BEGIN}, ${END}p" | while read sample
do
	INPUT_FASTQ_1=${INPUT_DIR}/'META20IMWJ27_'${sample}'_1.fastq.gz'
	INPUT_FASTQ_2=${INPUT_DIR}/'META20IMWJ27_'${sample}'_2.fastq.gz'
	OUTPUT_FASTQ_1=${OUTPUT_DIR}/${sample}'_fastp_1.fastq.gz'
	OUTPUT_FASTQ_2=${OUTPUT_DIR}/${sample}'_fastp_2.fastq.gz'
	OUTPUT_UNPAIRED=${OUTPUT_DIR}/${sample}'_unpaired.fastq.gz'
	fastp -i ${INPUT_FASTQ_1} -I ${INPUT_FASTQ_2} -o ${OUTPUT_FASTQ_1} -O ${OUTPUT_FASTQ_2} --thread 32 --unpaired1 ${OUTPUT_UNPAIRED} --unpaired2 ${OUTPUT_UNPAIRED} --correction --trim_poly_g --detect_adapter_for_pe
	#fastp -i ${INPUT_FASTQ_1} -I ${INPUT_FASTQ_2} -o ${OUTPUT_FASTQ_1} -O ${OUTPUT_FASTQ_2} --thread 32 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --unpaired1 ${OUTPUT_UNPAIRED} --unpaired2 ${OUTPUT_UNPAIRED}
done
conda deactivate
