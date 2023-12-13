# !/usr/bin/bash 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/21-2.spacers_blast_on_self_virus_plasmid.sh
# This script is for aliging recovered spacer sequences onto the viral and plasmid contigs predicted by geNomad, following script 7 and 23.
# Last modified: 23-10-3.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_QUERY_DIR=${BASE}/intermediates/7.crisprtools_extract
INPUT_REF_DIR=${BASE}/intermediates/23.genomad/total_pipeline
OUTPUT_DIR=${BASE}/intermediates/21.spacers_blast_on_targets/spacers_blast_on_mixed_virus_plasmid
mkdir -p ${OUTPUT_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi


# makeblastdb with mixed virus and plasmid contigs from all individuals.
total_virus_contigs=${OUTPUT_DIR}/blastdb/all_individuals/all_individuals_virus_contigs.fasta
total_plasmid_contigs=${OUTPUT_DIR}/blastdb/all_individuals/all_individuals_plasmid_contigs.fasta
mkdir -p ${OUTPUT_DIR}/blastdb/all_individuals
touch ${total_virus_contigs}
touch ${total_plasmid_contigs}
seq 1 10 | while read individual
do
        cat ${INPUT_REF_DIR}/${individual}/${individual}_noted.contigs_summary/${individual}_noted.contigs_virus.fna >> ${total_virus_contigs}
        cat ${INPUT_REF_DIR}/${individual}/${individual}_noted.contigs_summary/${individual}_noted.contigs_plasmid.fna >> ${total_plasmid_contigs}
done

        makeblastdb -in ${total_virus_contigs} -dbtype nucl -input_type fasta -out ${OUTPUT_DIR}/blastdb/all_individuals/all_individuals_virus_contigs -title all_individuals_virus_contigs
        makeblastdb -in ${total_plasmid_contigs} -dbtype nucl -input_type fasta -out ${OUTPUT_DIR}/blastdb/all_individuals/all_individuals_plasmid_contigs -title all_individuals_plasmid_contigs

seq ${BEGIN} ${END} | while read individual
do
	echo -----------------------------------------------------------------------------------
        echo ${individual}

	#blastn of spacers onto virus/plasmid contigs.
	cd ${INPUT_QUERY_DIR}
	ls ${individual}_*_*_spacers.fa | sed 's/_spacers.fa//g' | while read sample
	do
		QUERY=${INPUT_QUERY_DIR}/${sample}"_spacers.fa"
		
		REF_1=${OUTPUT_DIR}/blastdb/all_individuals/all_individuals_virus_contigs
		REF_2=${OUTPUT_DIR}/blastdb/all_individuals/all_individuals_plasmid_contigs

		OUTPUT_1=${OUTPUT_DIR}/${sample}"_all_individuals_virus_contigs.txt"
		OUTPUT_2=${OUTPUT_DIR}/${sample}"_all_individuals_plasmid_contigs.txt"

		echo ${QUERY}
		blastn -query ${QUERY}  -db ${REF_1}  -outfmt "6 qseqid sseqid pident nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus qseq sstrand frames " -num_threads 16 -max_target_seqs 5 -out ${OUTPUT_1}
		blastn -query ${QUERY}  -db ${REF_2}  -outfmt "6 qseqid sseqid pident nident qlen slen length mismatch positive ppos gapopen gaps qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus qseq sstrand frames " -num_threads 16 -max_target_seqs 5 -out ${OUTPUT_2}
	done
done
