
exome_dir <- "/Users/Eric/tcrseq/new/imseq/exome/"
rna_dir <- "/Users/Eric/tcrseq/new/imseq/rna/"
blood_dir <- "/Users/Eric/tcrseq/new/imseq/blood/"

file_exome <- list.files(path=exome_dir, pattern="*.tsv", full.names=T, recursive=FALSE)
file_rna <- list.files(path=rna_dir, pattern="*.tsv", full.names=T, recursive=FALSE)
file_blood <- list.files(path=blood_dir, pattern="*.tsv", full.names=T, recursive=FALSE)

results <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-13-2016.txt", sep="\t")
 
### get output of all clonotypes by patient id for blood, rna, and exome
clonotypes_all <- matrix(nrow=0, ncol=6)
clonotypes_exome <- matrix(nrow=0, ncol=6)
colnames(clonotypes_all) <- c("clonotype", "patient_id", "blood", "rna", "exome", "aa")
exome_clonotypes <- c()
rna_clonotypes <- c()
blood_clonotypes <- c()
for (i in 1:nrow(results))
{
  clonotypes_patient <- matrix(nrow=0, ncol=6)
  clonotypes_exome_patient <- matrix(nrow=0, ncol=6)
  exome_reads <- matrix(nrow=0, ncol=11)
  exome_file <- paste(results$exome_uuid[i], ".tsv", sep="")
  exome_path <- paste(exome_dir, exome_file, sep="")
  rna_reads <- matrix(nrow=0, ncol=11)
  rna_file <- paste(results$rna_uuid[i], ".tsv", sep="")
  rna_path <- paste(rna_dir, rna_file, sep="")
  blood_reads <- matrix(nrow=0, ncol=11)
  blood_file <- paste(results$blood_uuid[i], ".tsv", sep="")
  blood_path <- paste(blood_dir, blood_file, sep="")

  if(!file.exists(exome_path)) {next}
  exome_reads <- as.matrix(read.csv(exome_path, sep="\t", header=T))
  if(file.exists(rna_path)) {rna_reads <- as.matrix(read.csv(rna_path, sep="\t", header=T))
  } else {rna_reads <- NULL}
  if(file.exists(blood_path)) {blood_reads <- as.matrix(read.csv(blood_path, sep="\t", header=T))
  } else {blood_reads <- NULL}

  if (nrow(exome_reads)>0) 
  {
    for (j in 1:nrow(exome_reads))
    {
      cl <- exome_reads[j,10]
      aa <- exome_reads[j,11]
      match <- which(clonotypes_patient[,1]==cl)
      left_match <- exome_reads[j,4]
      right_match <- exome_reads[j,7]
      if (length(match)==0) 
      {
        clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$patient_uuid[i]), 0, 0, 1, aa))
#        clonotypes_exome_patient <- rbind(clonotypes_exome_patient, c(cl, aa, as.character(substr(results$sample_barcode[i], 1, 16)), 1, left_match, right_match))
      } else 
      {
        clonotypes_patient[match,5] <- as.character(as.numeric(clonotypes_patient[match,5]) + 1)
#        clonotypes_exome_patient[match,4] <- as.character(as.numeric(clonotypes_exome_patient[match,4]) + 1)
      }
      
    }    
  }
  
  if (length(rna_reads)>2)
  {
    for (j in 1:nrow(rna_reads))
    {
      cl <- rna_reads[j,10]
      aa <- rna_reads[j,11]
      match <- which(clonotypes_patient[,1]==cl)
      if (length(match)==0) {clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$patient_uuid[i]), 0, 1, 0, aa))  }
      else {clonotypes_patient[match,4] <- as.character(as.numeric(clonotypes_patient[match,4]) + 1)}
    }
  }

  if (length(blood_reads)>2)
  {
    for (j in 1:nrow(blood_reads))
    {
      cl <- blood_reads[j,10]
      aa <- blood_reads[j,11]
      match <- which(clonotypes_patient[,1]==cl)
      if (length(match)==0) {clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$patient_uuid[i]), 1, 0, 0, aa))  }
      else {clonotypes_patient[match,3] <- as.character(as.numeric(clonotypes_patient[match,3]) + 1)}
    }
  }

  if (is.null(rna_reads)) {clonotypes_patient[,4] <- rep(NA, nrow(clonotypes_patient))}
  if (is.null(blood_reads)) {clonotypes_patient[,3] <- rep(NA, nrow(clonotypes_patient))}
  clonotypes_all <- rbind(clonotypes_all, clonotypes_patient)
#  if (nrow(clonotypes_exome_patient)>0) {clonotypes_exome <- rbind(clonotypes_exome, clonotypes_exome_patient)}
}

write.table(clonotypes_all, "/Users/Eric/tcrseq/new/processed/pancan_clonotypes_all_7-14-2016.txt", sep="\t", quote=FALSE, row.names=FALSE)
