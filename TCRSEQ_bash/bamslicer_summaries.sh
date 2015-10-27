#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

base_dir="/scratch/data/bamslicer"
filename="${base_dir}/ref/BRCA.txt"
filelines=`cat $filename`
data_dir="${base_dir}/brca_slicer/data"
output_dir="${base_dir}/"
outdir="${base_dir}/data"
summary_dir="${base_dir}/summaries"
bam_outpath="/scratch"
NUM=0
QUEUE=""
MAX_NPROC=3

# declarations
region="1:1000000-1001000"
cghub_key=`cat /mnt/oncogxA/Administration/TCGA/cghub.key`

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

function bamslicer {
	URL="https://slicer.cghub.ucsc.edu/analyses"
	exome_id=$1
	ref="GRCh37-lite"
	format="bam"
	query_unmapped="${URL}/${exome_id}/slices?ref=${ref}&format=${format}&range=*"
	query_index="${URL}/${exome_id}/${exome_id}.bam.bai"
}

function fastq_convert {
	name=$1
	bamslicer $name
	curl -s "$query_region" -u ":${cghub_key}" > $bam_outpath/${name}.TCR.bam
	wait
	curl -s "$query_unmapped" -u ":${cghub_key}" > $bam_outpath/${name}.unmapped.bam
	wait
	local download_end=$(date +%s)
    local time_download=$(echo "$download_end - $time_start" | bc)
 	samtools index $bam_outpath/$name.TCR.bam
 	wait
 	samtools index $bam_outpath/$name.unmapped.bam
 	wait
 	samtools merge $bam_outpath/$name.TCRreg.bam $bam_outpath/$name.TCR.bam $bam_outpath/$name.unmapped.bam
 	wait
 	echo 'Overlap-Unmapped-count:' >> $bam_outpath/$name.summary.txt
 	samtools view -c $bam_outpath/$name.TCRreg.bam >> $bam_outpath/$name.summary.txt # unmapped-overlapping-count
 	wait
 	samtools sort $bam_outpath/$name.TCRreg.bam $bam_outpath/$name.TCRreg.sorted
	wait
	java -Xmx16g -jar /mnt/idash/Genomics/bin/picard-tools/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR INPUT=$bam_outpath/$name.TCRreg.sorted.bam OUTPUT=$bam_outpath/$name.TCRreg.sorted.dedup.bam METRICS_FILE=$bam_outpath/$name.metrics.txt
	wait
	bedtools bamtofastq -i $bam_outpath/$name.TCRreg.sorted.dedup.bam -fq $bam_outpath/$name.TCRreg.fastq
	wait
	bgzip $bam_outpath/$name.TCRreg.fastq
	wait
	local slice_end=$(date +%s)
	local time_slice=$(echo "$slice_end - $download_end" | bc)
	echo "$1	$time_download	$time_slice" >> ${output_dir}/summary.txt

	# move back
	mv $bam_outpath/$name.TCRreg.fastq.gz $outdir
	mv $bam_outpath/$name.summary.txt $outdir
#	mv $bam_outpath/$name.TCRreg.sorted.bam $outdir
	wait

	# clean up
	rm $bam_outpath/$name*
}

function main {
	name=$1
	bamslicer $name
	curl -s "$query_unmapped" -u ":${cghub_key}" > $bam_outpath/${name}.bam
	wait
	curl -s "$query_index" -u ":${cghub_key}" > $bam_outpath/${name}.bam.bai
	wait
	samtools idxstats $bam_outpath/${name}.bam > ${summary_dir}/${name}.txt
	wait
	rm $bam_outpath/$name*
}

# Main program
echo Start
while read f1 f2
do
	main $f1 $f2 &

    PID=$!
    queue $PID

    while [ $NUM -ge $MAX_NPROC ]; do
        checkqueue
        sleep 0.4
    done
done < $filename
wait # wait for all processes to finish before exit

