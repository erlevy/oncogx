#!/bin/bash
## Script to download and slice files
## Supports N processes at the same time
## For more processes change the max_jobs variable

base_dir="/scratch/data/bamslicer_exomes"
filename="${base_dir}/ref/BRCA.txt"
filelines=`cat $filename`
data_dir="${base_dir}/data"
output_dir="${base_dir}/"
outdir="${base_dir}/data"
summary_dir="${base_dir}/summaries"
bam_outpath="/scratch"
NUM=0
QUEUE=""
MAX_NPROC=10

# declarations
region_19="range=chr7:38295938-38407399&range=chr7:142000817-142510993&range=chr9:33618203-33662661&range=chr14:22090036-23014042"
region_37="range=7:38295938-38407399&range=7:142000817-142510993&range=9:33618203-33662661&range=14:22090036-23014042"
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
	ref=$2
	format="bam"
	query_region="${URL}/${exome_id}/slices?ref=${ref}&format=${format}&${region}"
	query_unmapped="${URL}/${exome_id}/slices?ref=${ref}&format=${format}&range=*"
	query_index="${URL}/${exome_id}/${exome_id}.bam.bai"
}

function fastq_convert {
	name=$1
	reference=$2
	bamslicer $name $reference
	local time_start=$(date +%s)
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
	samtools idxstats $bam_outpath/${name}.bam > ${summary_dir}/${name}.txt
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
#	local time_start=$(date +%s)
#	echo "Download start: " >> ${output_dir}/log.txt
#	echo $(date) >> ${output_dir}/log.txt
	fastq_convert $1 $2 &
#	wait
#	local download_end=$(date +%s)
#    local time_download=$(echo "$download_end - $time_start" | bc)
#	echo "Slicing start: " >> ${output_dir}/log.txt
#   echo $(date) >> ${output_dir}/log.txt
#	for file in $data_dir/$1/*.bam; do
#		[ ! -f "$file" ] && continue
#		echo $file
#		/home/edlevy/oncogx/Bam_to_Fastq.sh $file /scratch/data/Target-genesTR.ABDG-GenomicReg.bed ${output_dir}/slices /scratch/data/slices >> ${output_dir}/log.txt
#		wait
#	done
#	local slice_end=$(date +%s)
#	local time_slice=$(echo "$slice_end - $download_end" | bc)
#	echo "$1	$time_download	$time_slice" >> ${output_dir}/summary.txt
#	echo "Sample end: " >> ${output_dir}/log.txt
#	echo $(date) >> ${output_dir}/log.txt
#	rm $data_dir/$1/*.bam	# remove bam file for space after (didn't want to put an -r for the folder, figure out a safer way later)
}

# Main program
echo Start
echo "Exome	Download	Slice" > ${output_dir}/summary.txt
#echo "Log start: $(date)" > ${output_dir}/log.txt
while read f1 f2
do
#	echo $line >> ${output_dir}/log.txt
	main $f1 $f2 #>> ${output_dir}/log.txt & 

    PID=$!
    queue $PID

    while [ $NUM -ge $MAX_NPROC ]; do
        checkqueue
        sleep 0.4
    done
done < $filename
wait # wait for all processes to finish before exit
#echo "Log end: $(date)" >> ${output_dir}/log.txt
