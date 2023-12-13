# !/usr/bin/bash 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/21-1.spacers_blast_on_database.sh
# This script is for aliging recovered spacer sequences onto the virus reference constructed from public databases and viral contigs (vOTUs) from assembled metagenomic data, following script 7
# Last modified: 23-9-13.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/7.crisprtools_extract
OUTPUT_DIR=${BASE}/intermediates/21.spacers_blast_on_targets/spacers_blast_on_database
mkdir -p ${OUTPUT_DIR}

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
	echo -----------------------------------------------------------------------------------
	echo ${individual}
	cd ${INPUT_DIR}

	ls ${individual}_*_*_spacers.fa | sed 's/_spacers.fa//g' | while read sample
	do
		QUERY=${INPUT_DIR}/${sample}"_spacers.fa"
		
		REF_1=~/crisprome/virome_database/blastdb/HuVirDB
		REF_2=~/crisprome/virome_database/blastdb/virus_sequences
		REF_3=~/crisprome/virome_database/blastdb/plasmid_seq

		OUTPUT_1=${OUTPUT_DIR}/${sample}"_HuVirDB.txt"
		OUTPUT_2=${OUTPUT_DIR}/${sample}"_virus_sequences.txt"
		OUTPUT_3=${OUTPUT_DIR}/${sample}"_plasmid_seq.txt"

		echo ${QUERY}
		blastn -query ${QUERY}  -db ${REF_1}  -outfmt "6 qseqid sseqid pident nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus qseq sstrand frames " -num_threads 16 -max_target_seqs 5 -out ${OUTPUT_1}
		blastn -query ${QUERY}  -db ${REF_2}  -outfmt "6 qseqid sseqid pident nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus qseq sstrand frames " -num_threads 16 -max_target_seqs 5 -out ${OUTPUT_2}
		blastn -query ${QUERY}  -db ${REF_3}  -outfmt "6 qseqid sseqid pident nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus qseq sstrand frames " -num_threads 16 -max_target_seqs 5 -out ${OUTPUT_3}
	done
done
