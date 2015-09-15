#!/bin/bash


name=`basename $1 ".fastq.gz"`
dir=`dirname $1`

/scratch/clonotyper_hs/scripts/clonotypeR detect $dir/$name.fastq.gz

wait

/scratch/clonotyper_hs/scripts/clonotypeR extract $dir/$name.fastq.gz 

wait




	
