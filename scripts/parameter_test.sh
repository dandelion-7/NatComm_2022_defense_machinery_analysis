# !/usr/bin/bash

if [ -z "$1" ]; then
	echo "there's no input parameter"
else
	for parameter in $*
       	do
		echo the input parameter is ${parameter}
	done
fi

INPUT_ILLUMINA_DIR=~/public_raw_data/wangjun_2022_natComm_ngs+tgs_SV/data/illumina/time_series
if [ -z "$1" ]; then
	ls ${INPUT_ILLUMINA_DIR}
else
	ls ${INPUT_ILLUMINA_DIR} | sed -n "${1}, ${2}p"
fi
