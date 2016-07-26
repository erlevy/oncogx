input='/scratch/immunoseq/fastq'
ref='/mnt/idash/Genomics/data_resources/references-and-indexes/hg19/bwa_v7'
output='/scratch/immunoseq/aligned'
exome='/mnt/idash/Genomics/data_resources/exome/baits_targets'

/mnt/idash/Genomics/bin/bwa mem -M -t 26 ${ref}/hg19_lite.fa ${input}/merged_R1.fastq ${input}/merged_R2.fastq > ${output}/aln.sam
wait

samtools view -Sb ${output}/aln.sam > ${output}/aln.bam
wait

samtools sort -@ 28 ${output}/aln.bam ${output}/aln.sorted
wait

java -Xmx16g -jar /mnt/idash/Genomics/bin/picard-tools/MarkDuplicates.jar INPUT=${output}/aln.sorted.bam OUTPUT=${output}/aln.sorted.dedup.bam METRICS_FILE=${output}/metrics.txt REMOVE_DUPLICATES=TRUE

java -Xmx16g -jar /mnt/idash/Genomics/bin/picard-tools/CalculateHsMetrics.jar  BAIT_INTERVALS=${exome}/humanV4-baits_hg19_lite.interval_list TARGET_INTERVALS=${exome}/humanV4-targets_hg19_lite.interval_list INPUT=${output}/aln.sorted.dedup.bam OUTPUT=${output}/aln.sorted.dedup.bam.metrics

samtools index ${output}/aln.sorted.bam
wait

