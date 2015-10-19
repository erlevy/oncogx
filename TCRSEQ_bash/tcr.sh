### mitcr ###

input_bam='/scratch/immunoseq/sliced/dedup/aln.sorted.dedup.TCRreg.fastq.gz'
output_mitcr='/scratch/immunoseq/results/mitcr'

# mitcr <options> <input file name> <output file name>
mitcr -pset flex $input_bam ${output_mitcr}/result.txt

### imseq ###

imseq_path='/home/edlevy/imseq_1.0.1-linux64'
output_imseq='/scratch/immunoseq/results/imseq'

# mitcr <options> <input file name> <output file name>
${imseq_path}/imseq -ref ${imseq_path}/segment-reference.fa -o ${output_imseq}/imseq_output.tsv $input_bam