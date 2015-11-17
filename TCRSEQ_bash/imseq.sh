#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

# need to run program from within imseq directory
#imseq -ref segment-reference.fa -o output-file.tsv input-file.fastq.gz

workdir=$1
CNT=1
outdir="/scratch/data/bamslicer_exomes/imseq"
for file in $workdir/*fastq; do

  name=`basename $file ".TCRreg.fastq"`
  ./imseq -ref Homo.Sapiens.TRB.fa -o $outdir/$name.tsv $file &
  
  if [ $(($CNT % 26)) -eq 0 ]; then
    echo "##DOING 26##"
    wait
  fi
  
  let CNT+=1

done