#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

filename='/scratch/data/brca_slicer/ref/filtered_exomes_brca_25.txt'
filelines=`cat $filename`
data_dir='/scratch/data/brca_slicer/data'
output_dir='/scratch/data/brca_slicer/'
NUM=0
QUEUE=""
MAX_NPROC=10

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}

function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
        if [ -d /proc/$PID  ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}

function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
        if [ ! -d /proc/$PID ] ; then
            regeneratequeue # at least one PID has finished
            break
        fi
    done
}

function main {
	local time_start=$(date +%s)
	echo "Download start: " >> ${output_dir}/log.txt
	echo $(date) >> ${output_dir}/log.txt
    gtdownload -c /mnt/oncogxA/Administration/TCGA/cghub.key --max-children 3 -p ${output_dir}/data -d $1 >> ${output_dir}/log.txt
    wait
	local download_end=$(date +%s)
    local time_download=$(echo "$download_end - $time_start" | bc)
	echo "Slicing start: " >> ${output_dir}/log.txt
    echo $(date) >> ${output_dir}/log.txt
	for file in $data_dir/$1/*.bam; do
		[ ! -f "$file" ] && continue
		echo $file
		/home/edlevy/oncogx/Bam_to_Fastq.sh $file /scratch/data/Target-genesTR.ABDG-GenomicReg.bed ${output_dir}/slices /scratch/data/slices >> ${output_dir}/log.txt
		wait
	done
	local slice_end=$(date +%s)
    local time_slice=$(echo "$slice_end - $download_end" | bc)
	echo "$1	$time_download	$time_slice" >> ${output_dir}/summary.txt
	echo "Sample end: " >> ${output_dir}/log.txt
	echo $(date) >> ${output_dir}/log.txt
	rm $data_dir/$1/*.bam	# remove bam file for space after (didn't want to put an -r for the folder, figure out a safer way later)
}

# Main program
echo Start
echo "Exome	Download	Slice" > ${output_dir}/summary.txt
echo "Log start: $(date)" > ${output_dir}/log.txt
for line in $filelines ; do
	echo $line >> ${output_dir}/log.txt
	main $line >> ${output_dir}/log.txt & 

    PID=$!
    queue $PID

    while [ $NUM -ge $MAX_NPROC ]; do
        checkqueue
        sleep 0.4
    done
done
wait # wait for all processes to finish before exit
echo "Log end: $(date)" >> ${output_dir}/log.txt
