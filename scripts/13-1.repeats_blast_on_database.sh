# !/usr/bin/bash
#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/13-1.repeats_blast_on_database.sh

# This script is for blast all the recovered repeats against the blast database constructed with repeats sequences from five studies, following script 7
# Last modified: 23-7-13.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/7.crisprtools_extract
OUTPUT_DIR=${BASE}/intermediates/13.repeats_blast_on_database/blastn_results
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
	echo ${individual}
	ls ${INPUT_DIR} | grep "repeats.fa" | grep "${individual}_.*_TD" | sed "s/.fa//g"  | while read sample
do
	query=${INPUT_DIR}/${sample}".fa"
	echo ${query}
	db=~/crisprome/crispr_cas_database/summary/direct_repeats/blastdb/dedup_dr_summary.fasta
	OUTPUT_FILE=${OUTPUT_DIR}/${sample}"_blastn_database.txt"	

	blastn -query ${query}  -db ${db}  -outfmt "6 qseqid sseqid pident nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus qseq sstrand frames " -num_threads 10 -max_target_seqs 1 -out ${OUTPUT_FILE}
done
done

echo "qseqid	sseqid	pident	nident	qlen	slen	length	mismatch	positive	ppos	gapopen	gaps	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp	qcovus	qseq	sstrand	frames" > ${OUTPUT_DIR}/total_blastn_database.txt
cat ${OUTPUT_DIR}/*TD*.txt >> ${OUTPUT_DIR}/total_blastn_database.txt
