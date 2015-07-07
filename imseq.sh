#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

# need to run program from within imseq directory
imseq -ref segment-reference.fa -o output-file.tsv input-file.fastq.gz