#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

# program is in path
#mitcr ­pset flex in.fastq.gz result.txt

workdir=$1
CNT=1
outdir="/scratch/data/bamslicer/mitcr"
for file in $workdir/*fastq.gz; do

  name=`basename $file ".TCRreg.fastq.gz"`
  /usr/bin/mitcr ­pset flex $file $outdir/$name.txt &

  
  if [ $(($CNT % 26)) -eq 0 ]; then
    echo "##DOING 26##"
    wait
  fi
  
  let CNT+=1
done
