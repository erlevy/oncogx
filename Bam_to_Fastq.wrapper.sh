#!/bin/bash

bam_dir='/scratch/immunoseq/aligned'
out_dir='/scratch/immunoseq/sliced/dedup'
CNT=1

for file in ${bam_dir}/*.bam; do
	name=`basename $file ".bam"`
	dir=`dirname $file`
	
	/home/edlevy/oncogx/Bam_to_Fastq.sh $file /scratch/data/Target-genesTR.ABDG-GenomicReg.bed ${out_dir} ${out_dir} & 

	if [ $(($CNT % 26)) -eq 0 ]; then
		echo "##DOING 26##"
		wait
	fi
	let CNT+=1
done
