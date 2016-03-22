
imseq_dir <- '/Users/Eric/tcga/clonotypes/final/blood/'

results <- read.csv("blood_normal_results_with_clonotypes.csv")
exome_freq <- as.numeric(results$normal_cdr3)/as.numeric(results$normal_reads)
results <- cbind(results, exome_freq)

# comparisons
results2 <- read.csv("cdr3_results_with_groups_complete.txt", sep="\t")
remove <- c()
exome_count <- c()
for (i in 1:nrow(results))
{
  patient <- as.character(results$patient_id[i])
  matching <- which(as.character(results2$patient_id)==patient)
  if (length(matching)==0) {remove <- c(remove, i)}
  else {exome_count <- c(exome_count, results2$exome_reads[matching])}
}

# need to make a final set that has blood
blood_cdr3 <- c()
blood_reads <- c()
blood_id <- c()
matching_reads <- c()
blood_clonotypes <- c()
for (i in 1:nrow(results2))
{
  patient <- as.character(results2$patient_id[i])
  matching <- which(as.character(results$patient_id)==patient)
  if (length(matching)==0) 
  {
    blood_cdr3 <- c(blood_cdr3, NA)
    blood_reads <- c(blood_reads, NA)
    blood_id <- c(blood_id, NA)
    matching_reads <- c(matching_reads, NA)
    blood_clonotypes <- c(blood_clonotypes, NA)
  } else 
  {
    blood_cdr3 <- c(blood_cdr3, results$normal_cdr3[matching])
    blood_reads <- c(blood_reads, results$normal_reads[matching])
    blood_id <- c(blood_id, as.character(results$normal_id[matching]))
    matching_reads <- c(matching_reads, as.character(results$matching[matching]))
    blood_clonotypes <- c(blood_clonotypes, as.character(results$normal_clonotypes[matching]))
  }  
}
results3 <- cbind(results2, blood_id, blood_reads, blood_cdr3, blood_clonotypes, matching_reads)

#write.table(results3, "cdr3_results_with_blood_1078.txt", quote=FALSE, sep="\t", row.names=FALSE)

# calculate stuff
results <- results3

exome_rpm <- as.numeric(results$exome_imseq)/as.numeric(results$exome_reads)*1000000
rna_rpm <- as.numeric(results$rna_imseq)/as.numeric(results$rna_reads)*1000000
blood_rpm <- as.numeric(results$blood_cdr3)/as.numeric(results$blood_reads)*1000000

have_match <- which(as.numeric(as.character(results$matching_reads))>0)
no_match <- which(as.numeric(as.character(results$matching_reads))==0)

# get shared clonotype in the blood samples
# imseq: col 1: read ID, col 4: V match, col 7: J match, col 10: nucleotide, col 11: AA
shared <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
imseq <- matrix(nrow=0, ncol=3)
colnames(imseq) <- c("nucleotide_sequence", "aa_sequence", "number_of_patients")
blood_shared <- rep(0, nrow(results))
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  imseq_file <- paste(results$blood_id[i], ".tsv", sep="")
  imseq_path <- paste(imseq_dir, imseq_file, sep="")
  if (!file.exists(imseq_path))
  {
    next
  }
  if(file.info(imseq_path)$size!=0)
  {
    imseq_reads <- as.matrix(read.csv(imseq_path, sep="\t", header=T))
    imseq_reads <- imseq_reads[!duplicated(imseq_reads[,10]),, drop=FALSE] 
  }
  if (nrow(imseq_reads)>0)
  {
    for (j in 1:nrow(imseq_reads))
    {
      nucleotide <- imseq_reads[j,10]
      if (nucleotide==shared) {blood_shared[i] <- 1}
      aa <- imseq_reads[j,11]
      imseq_match <- which(imseq[,1]==nucleotide)
      if (length(imseq_match)==0)
      {
        imseq <- rbind(imseq, c(nucleotide, aa, "1"))
      } else
      {
        imseq[imseq_match,3] <- as.character(as.numeric(imseq[imseq_match,3])+1)
      }
    }
  }
}
imseq_blood <- imseq
write.table(imseq, "repeated_clonotypes_normal.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# get shared clonotype in the tumor samples
shared <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
imseq <- matrix(nrow=0, ncol=3)
colnames(imseq) <- c("nucleotide_sequence", "aa_sequence", "number_of_patients")
tumor_shared <- rep(0, nrow(results))
imseq_dir <- '/Users/Eric/tcga/clonotypes/final/exome/'
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  imseq_file <- paste(results$exome_id[i], ".tsv", sep="")
  imseq_path <- paste(imseq_dir, imseq_file, sep="")
  if (!file.exists(imseq_path))
  {
    next
  }
  if(file.info(imseq_path)$size!=0)
  {
    imseq_reads <- as.matrix(read.csv(imseq_path, sep="\t", header=T))
    imseq_reads <- imseq_reads[!duplicated(imseq_reads[,10]),, drop=FALSE] 
  }
  if (nrow(imseq_reads)>0)
  {
    for (j in 1:nrow(imseq_reads))
    {
      nucleotide <- imseq_reads[j,10]
      if (nucleotide==shared) {tumor_shared[i] <- 1}
      aa <- imseq_reads[j,11]
      imseq_match <- which(imseq[,1]==nucleotide)
      if (length(imseq_match)==0)
      {
        imseq <- rbind(imseq, c(nucleotide, aa, "1"))
      } else
      {
        imseq[imseq_match,3] <- as.character(as.numeric(imseq[imseq_match,3])+1)
      }
    }
  }
}
imseq_tumor <- imseq
write.table(imseq, "repeated_clonotypes_tumor.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

results3 <- cbind(results3, tumor_shared, blood_shared)

length(which(as.numeric(as.character(results3$matching_reads))>0)) # 34
length(which(as.numeric(as.character(results3$tumor_shared))>0)) # 62
length(which(as.numeric(as.character(results3$blood_shared))>0)) # 72
length(intersect(which(as.numeric(as.character(results3$blood_shared))>0), which(as.numeric(as.character(results3$tumor_shared))>0))) # 33
