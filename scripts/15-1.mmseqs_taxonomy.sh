# !usr/bin/bash
#########################################################################################################################################################################
# Script: ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/15.mmseqs_taxonomy.sh
# This script is for assigning taxonomy of megahit-assembled contigs with mmseqs2, following script 10.
# Last modified: 23-8-28.
#########################################################################################################################################################################

BASE=~/crisprome/human_temporal_metagenome_wangjun_2022natComm
INPUT_DIR=${BASE}/intermediates/10.megahit
OUTPUT_TAXONOMY_DIR=${BASE}/intermediates/15.mmseqs_taxonomy/taxonomy_new
OUTPUT_QUERYDB_DIR=${BASE}/intermediates/15.mmseqs_taxonomy/querydb_new
mkdir -p ${OUTPUT_TAXONOMY_DIR}
mkdir -p ${OUTPUT_QUERYDB_DIR}

# Divide all the files into multiple parallel tasks.
if [ -z "$1" ]; then
	BEGIN=1
	END=10
else
	BEGIN=${1}
	END=${2}
fi

source activate MMseqs2
seq ${BEGIN} ${END} | while read i
do
	cd ${OUTPUT_QUERYDB_DIR}
	FASTA=${INPUT_DIR}/${i}/${i}".contigs.fa"
	FASTA_NOTED=${INPUT_DIR}/${i}/${i}"_noted.contigs.fa"
	sed -e "s/>/>${i}_/g" -e "s/ /_/g" ${FASTA} > ${FASTA_NOTED}
	echo ${FASTA_NOTED}
	mmseqs createdb ${FASTA_NOTED} ${i} # querydb is removed after assigning taxonomies to save storage.

	QUERYDB=${OUTPUT_QUERYDB_DIR}/${i}
	UNIREFDB=~/genome/uniref/mmseqs_uniref100/mmseqs_uniref100_db
	tmp=${OUTPUT_TAXONOMY_DIR}/${i}_tmp
	mkdir -p ${tmp}
	RESULT=${OUTPUT_TAXONOMY_DIR}/${i}
	TSV=${OUTPUT_TAXONOMY_DIR}/${i}.txt

	mmseqs taxonomy ${QUERYDB} ${UNIREFDB} ${RESULT} ${tmp} --lca-ranks species,genus,family,order,class,phylum,kingdom
	mmseqs createtsv ${QUERYDB} ${RESULT} ${TSV}
done
conda deactivate
