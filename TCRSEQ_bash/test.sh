#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

base_dir="/scratch/data/bamslicer"
filename="${base_dir}/ref/BRCA_2.txt"
filelines=`cat $filename`
data_dir="${base_dir}/brca_slicer/data"
output_dir="${base_dir}/"
outdir="${base_dir}/data"
bam_outpath="/scratch"
cghub_key=`cat /mnt/oncogxA/Administration/TCGA/cghub.key`
NUM=0
QUEUE=""
MAX_NPROC=10

# declarations

function bamslicer {
	URL="https://slicer.cghub.ucsc.edu/analyses"
	exome_id="$1"
	ref=$2
	format="bam"
	check_reference $1 $2
	query_region="${URL}/${exome_id}/slices?ref=${ref}&format=${format}&${region}"
	query_unmapped="${URL}/${exome_id}/slices?ref=${ref}&format=${format}&range=*"
}

function check_reference {
	query_header="https://slicer.cghub.ucsc.edu/analyses/${1}/header?ref=${2}"
	region="range=7:38295938-38407399&range=7:142000817-142510993&range=9:33618203-33662661&range=14:22090036-23014042"
	region="range=chr7:38295938-38407399&range=chr7:142000817-142510993&range=chr9:33618203-33662661&range=chr14:22090036-23014042"
}

function main {
	bamslicer $1 $2
	echo "$query_header"
	curl -s "$query_header" -u ":${cghub_key}" > $bam_outpath/${1}.header.bam

}

while read f1 f2
do
	main $f1 $f2
done < $filename