#!/bin/bash

bam_dir='/scratch/data/brca_slicer/slices'
out_dir='/scratch/data/brca_slicer/fastq'
CNT=1

function main {
	bedtools bamtofastq -i $file -fq $out_dir/$name.TCRreg.fastq
	wait
	bgzip $out_dir/$name.TCRreg.fastq
	wait
}


for file in ${bam_dir}/*.bam; do
	name=`basename $file ".bam"`
	dir=`dirname $file`
	#/home/edlevy/scripts/Bam_to_Fastq-fromAnnai-Single.sh $file /scratch/data/Target-genesTR.ABDG-GenomicReg.bed /scratch/data/test/100/fastq /scratch/data/test/100/tcr_bam & 

	main $file $name $dir &

	if [ $(($CNT % 26)) -eq 0 ]; then
		echo "##DOING 26##"
		wait
	fi
	let CNT+=1
done
