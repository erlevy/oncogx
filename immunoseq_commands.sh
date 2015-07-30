
mitcr -pset flex -quality 0 /scratch/immunoseq/sliced/dedup/aln.sorted.dedup.TCRreg.fastq.gz /scratch/immunoseq/results/mitcr/mitcr_results_quality_0.txt

./imseq -ref Homo.Sapiens.TRB.fa -rlog /scratch/immunoseq/results/imseq/imseq_rejected.tsv -o /scratch/immunoseq/results/imseq/imseq_single.tsv /scratch/immunoseq/sliced/dedup/aln.sorted.dedup.TCRreg.fastq
