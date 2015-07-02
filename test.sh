#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

base_dir="/scratch/data/bamslicer"
filename="${base_dir}/ref/BRCA.txt"
filelines=`cat $filename`
data_dir="${base_dir}/brca_slicer/data"
output_dir="${base_dir}/"
outdir="${base_dir}/data"
bam_outpath="/scratch"
NUM=0
QUEUE=""
MAX_NPROC=10

# declarations
region="range=7:38295938-38407399&range=7:142000817-142510993&range=9:33618203-33662661&range=14:22090036-23014042"

while read f1 f2
do
	echo $f1
	echo $f2
done < $filename