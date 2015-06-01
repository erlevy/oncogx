input='/scratch/immunoseq/fastq'
ref='/mnt/idash/Genomics/data_resources/references-and-indexes/hg19/bwa_v7'
output='/scratch/immunoseq/aligned'

/mnt/idash/Genomics/bin/bwa mem -M -t 26 ${ref}/hg19_lite.fa ${input}/merged_R1.fastq ${input}/merged_R2.fastq > ${output}/aln.sam