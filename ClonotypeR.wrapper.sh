#!/bin/bash

workdir=$1
CNT=1
scriptsDir="/home/edlevy/oncogx"
#source ~/.bashrc
for file in $workdir/*fastq.gz; do

  echo $file
  ${scriptsDir}/ClonotypeR.sh  $file &
  

  
  if [ $(($CNT % 26)) -eq 0 ]; then
    echo "##DOING 26##"
    wait
  fi
  
  let CNT+=1
done
