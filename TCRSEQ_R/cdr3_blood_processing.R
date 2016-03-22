
normal_dir <- '/Users/Eric/tcga/clonotypes/final/blood/'
summaries_normal <- '/Users/Eric/tcga/clonotypes/final/blood_summaries/'
exome_dir <- '/Users/Eric/tcga/clonotypes/final/exome/'
summaries_exome <- '/Users/Eric/tcga/clonotypes/final/exome_summaries/'
#clinical_dir <- 'clinical/'
reference_dir <- '/Users/Eric/tcga/results/'

#######
# IMSEQ Results Files
#######

file_summaries_normal <- list.files(path=summaries_normal, pattern="*.txt", full.names=T, recursive=FALSE)
file_normal <- list.files(path=normal_dir, pattern="*.tsv", full.names=T, recursive=FALSE)

file_summaries_exome <- list.files(path=summaries_exome, pattern="*.txt", full.names=T, recursive=FALSE)
file_exome <- list.files(path=exome_dir, pattern="*.tsv", full.names=T, recursive=FALSE)

#lymphocytes <- read.csv(paste(clinical_dir, 'lymphocyte_data_table.csv', sep=""), header=TRUE)
#clinical <- read.csv(paste(clinical_dir, 'clinical_brca.txt', sep=""), header=TRUE, sep="\t")
normal_samples <- read.csv(paste(reference_dir, 'lymphocyte_data_table_normal.csv', sep=""), header=TRUE)
tumor_samples <- read.csv(paste(reference_dir, 'lymphocyte_data_table.csv', sep=""), header=TRUE)

results <- matrix(nrow=length(file_summaries_exome), ncol=8)
cols <- c("patient_id", "tumor_id", "normal_id", "normal_reads", "normal_cdr3", "exome_cdr3", "matching", "normal_clonotypes")# , "days", "vital")
colnames(results) <- cols
# data matrix: patient id, sample id, exome id, exome total read count, exome imseq read count, rna id, rna total read count,
#              rna imseq read count, lymphocyte percentage, presence of lymphocytes, days to death/last contact, vital status

# put in tumor id
for (i in 1:length(file_summaries_exome))
{
  summary <- read.table(file_summaries_exome[i])
  base <- basename(file_summaries_exome[i])
  id <- substring(base, first=1, last=nchar(base)-4)
  results[i,2] <- id
}

# put in patient, normal ids, and normal read count
for (i in 1:nrow(results))
{
  tumor <- results[i,2]
  tumor_row <- which(tumor_samples[,9]==tumor)[1]
  patient <- as.character(tumor_samples[tumor_row, 1])
  results[i,1] <- patient
  normal_row <- which(normal_samples[,1]==patient)
  if (length(normal_row)!=0) 
  {
    normal_row <- normal_row[1]
    normal_id <- as.character(normal_samples[normal_row, 9])
    normal_file <- paste(summaries_normal, normal_id, sep="")
    normal_file <- paste(normal_file, ".txt", sep="")
    summary <- read.table(normal_file)
    reads <- sum(c(summary[,3], summary[,4]))
  } else 
  {
    normal_id <- NA
    reads <- 0
  }
  results[i,3] <- normal_id
  results[i,4] <- reads  
}

# put in imseq read counts for normal
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  normal_clonotypes <- 0
  imseq_file <- paste(results[i,3], ".tsv", sep="")
  imseq_path <- paste(normal_dir, imseq_file, sep="")
  if(file.exists(imseq_path))
  {
    if(file.info(imseq_path)$size!=0)
    {
      imseq_reads_file <- read.csv(imseq_path, sep="\t", header=T)
      imseq_reads <- nrow(imseq_reads_file)
      imseq_reads_file <- as.matrix(imseq_reads_file)
      imseq_reads_file <- imseq_reads_file[!duplicated(imseq_reads_file[,10]),, drop=FALSE] 
      normal_clonotypes <- nrow(imseq_reads_file)
    } else
    {
      imseq_reads <- NA
      normal_clonotypes <- NA
    }
  } else
  {
    imseq_reads <- NA
    normal_clonotypes <- NA
  }
  results[i,5] <- imseq_reads
  results[i,8] <- normal_clonotypes
}

# put in imseq read counts for tumor
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  imseq_file <- paste(results[i,2], ".tsv", sep="")
  imseq_path <- paste(exome_dir, imseq_file, sep="")
  if(file.exists(imseq_path))
  {
    if(file.info(imseq_path)$size!=0)
    {
      imseq_reads <- nrow(read.csv(imseq_path, sep="\t", header=T))
    } else
    {
      imseq_reads <- NA
    }
  } else
  {
    imseq_reads <- NA
  }
  results[i,6] <- imseq_reads
}

# get matching
for (i in 1:nrow(results))
{
  if (is.na(results[i,3]))
  {
    num_matching <- NA
  } else
  {
    normal_file <- paste(results[i,3], ".tsv", sep="")
    normal_path <- paste(normal_dir, normal_file, sep="")
    tumor_file <- paste(results[i,2], ".tsv", sep="")
    tumor_path <- paste(exome_dir, tumor_file, sep="")
    normal_reads <- read.csv(normal_path, sep="\t", header=T)
    tumor_reads <- read.csv(tumor_path, sep="\t", header=T)
    num_matching <- 0
    if (nrow(tumor_reads)>0)
    {
      for (j in 1:nrow(tumor_reads))
      {
        nuc <- as.character(tumor_reads$cdrNucSeq[j])
        normal_match <- which(normal_reads$cdrNucSeq==nuc)
        if (length(normal_match)>0) {num_matching <- num_matching + 1}
      }
    }    
  }
  results[i,7] <- num_matching
} 
  
# remove matching since i dont have it yet
results <- na.omit(results)
#results <- results[,-6]

write.table(results, "blood_normal_results_with_clonotypes.csv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")

# check shared
repeated <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
patients_normal <- c()
patients_tumor <- c()
for (i in 1:nrow(results))
{
  normal_file <- paste(results[i,3], ".tsv", sep="")
  normal_path <- paste(normal_dir, normal_file, sep="")
  imseq <- read.csv(normal_path, sep="\t", header=T)
  matching <- which(imseq$cdrNucSeq==repeated)
  if (length(matching)>0)
  {
    patients_normal <- c(patients_normal, results[i,1])
  }
  tumor_file <- paste(results[i,2], ".tsv", sep="")
  tumor_path <- paste(exome_dir, tumor_file, sep="")
  imseq <- read.csv(tumor_path, sep="\t", header=T)
  matching <- which(imseq$cdrNucSeq==repeated)
  if (length(matching)>0)
  {
    patients_tumor <- c(patients_tumor, results[i,1])
  }
}
intersection <- intersect(patients_tumor, patients_normal)

# counting
exome <- as.numeric(results[,6])
normal <- as.numeric(results[,5])
total <- sum(exome>0 & normal>0)
