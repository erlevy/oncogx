#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

# need to run program from within imseq directory
#imseq -ref segment-reference.fa -o output-file.tsv input-file.fastq.gz

workdir=$1
CNT=1
outdir="/scratch/data/bamslicer/imseq"
datadir="/scratch/data/bamslicer/data"
for file in $workdir/*fastq.gz; do

  name=`basename $file ".TCRreg.fastq.gz"`
  gzip -d $file
  wait
  ./imseq -ref Homo.Sapiens.TRB.fa -o $outdir/$name.tsv $datadir/$name.TCRreg.fastq &
  
  if [ $(($CNT % 26)) -eq 0 ]; then
    echo "##DOING 26##"
    wait
  fi
  
  let CNT+=1
  
  gzip $file
  wait
done