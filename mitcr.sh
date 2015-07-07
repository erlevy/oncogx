#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

# program is in path
mitcr ­pset flex in.fastq.gz result.txt
