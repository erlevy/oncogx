#!/bin/bash


name=`basename $1 ".fastq.gz"`
dir=`dirname $1`

/home/edlevy/clonotyper_hs/scripts/clonotypeR detect $dir/$name.fastq.gz

wait

/home/edlevy/clonotyper_hs/scripts/clonotypeR extract $dir/$name.fastq.gz 

wait




	
