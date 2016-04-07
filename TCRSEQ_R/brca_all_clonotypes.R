
exome_dir <- '/Users/Eric/tcga/clonotypes/final/exome/'
rna_dir <- '/Users/Eric/tcga/clonotypes/final/rna/'
blood_dir <- '/Users/Eric/tcga/clonotypes/final/blood/'

#######
# IMSEQ Results Files
#######

file_exome <- list.files(path=exome_dir, pattern="*.tsv", full.names=T, recursive=FALSE)
file_rna <- list.files(path=rna_dir, pattern="*.tsv", full.names=T, recursive=FALSE)
file_blood <- list.files(path=blood_dir, pattern="*.tsv", full.names=T, recursive=FALSE)

results <- read.csv("/Users/Eric/cdr3_results_with_expression_blood_ptprc_1078.txt", sep="\t")

clonotypes_all <- matrix(nrow=0, ncol=6)
clonotypes_exome <- matrix(nrow=0, ncol=6)
colnames(clonotypes_all) <- c("clonotype", "patient_id", "blood", "rna", "exome")
exome_clonotypes <- c()
rna_clonotypes <- c()
blood_clonotypes <- c()
for (i in 1:nrow(results))
{
  clonotypes_patient <- matrix(nrow=0, ncol=6)
  clonotypes_exome_patient <- matrix(nrow=0, ncol=6)
  exome_reads <- matrix(nrow=0, ncol=11)
  exome_file <- paste(results$exome_id[i], ".tsv", sep="")
  exome_path <- paste(exome_dir, exome_file, sep="")
  rna_reads <- matrix(nrow=0, ncol=11)
  rna_file <- paste(results$rna_id[i], ".tsv", sep="")
  rna_path <- paste(rna_dir, rna_file, sep="")
  blood_reads <- matrix(nrow=0, ncol=11)
  blood_file <- paste(results$blood_id[i], ".tsv", sep="")
  blood_path <- paste(blood_dir, blood_file, sep="")
  
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
        clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$patient_id[i]), 0, 0, 1, aa))
        clonotypes_exome_patient <- rbind(clonotypes_exome_patient, c(cl, aa, as.character(substr(results$sample_id[i], 1, 16)), 1, left_match, right_match))
      } else 
      {
        clonotypes_patient[match,5] <- as.character(as.numeric(clonotypes_patient[match,5]) + 1)
        clonotypes_exome_patient[match,4] <- as.character(as.numeric(clonotypes_exome_patient[match,4]) + 1)
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
      if (length(match)==0) {clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$patient_id[i]), 0, 1, 0, aa))  }
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
      if (length(match)==0) {clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$patient_id[i]), 1, 0, 0, aa))  }
      else {clonotypes_patient[match,3] <- as.character(as.numeric(clonotypes_patient[match,3]) + 1)}
    }    
  }
  
  if (is.null(rna_reads)) {clonotypes_patient[,4] <- rep(NA, nrow(clonotypes_patient))}
  if (is.null(blood_reads)) {clonotypes_patient[,3] <- rep(NA, nrow(clonotypes_patient))}
  clonotypes_all <- rbind(clonotypes_all, clonotypes_patient)
  clonotypes_exome <- rbind(clonotypes_exome, clonotypes_exome_patient)
}

write.table(clonotypes_all, "clonotypes_all_with_aa.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(clonotypes_exome, "clonotypes_exome.txt", sep="\t", quote=FALSE, row.names=FALSE)

clonotypes_small <- clonotypes_all[,c(2,4)]
clonotypes_small <- clonotypes_small[!duplicated(clonotypes_small[,1]),]
setdiff(results$patient_id[which(is.na(results$blood_id))], clonotypes_small[which(is.na(clonotypes_small[,2])),1])

normal_clonotypes <- read.csv("repeated_clonotypes_normal.txt", sep="\t")
tumor_clonotypes <- read.csv("repeated_clonotypes_tumor.txt", sep="\t")

clonotypes_small <- clonotypes_all[which(as.numeric(clonotypes_all[,5])>0),1]
clonotypes_small <- sort(unique(clonotypes_small))
  
results_complete <- cbind(results[,1:5], exome_clonotypes, results[,6:8], rna_clonotypes, results[,9:12])
cols <- c("patient_id", "sample_id", "exome_id", "exome_reads", "exome_imseq", "exome_clonotypes", "rna_id", "rna_reads", 
          "rna_imseq", "rna_clonotypes", "lymphocyte_percent", "lymphocyte_presence", "days", "vital")
colnames(results_complete) <- cols

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
  }
  else 
  {
    blood_cdr3 <- c(blood_cdr3, results$normal_cdr3[matching])
    blood_reads <- c(blood_reads, results$normal_reads[matching])
    blood_id <- c(blood_id, as.character(results$normal_id[matching]))
    matching_reads <- c(matching_reads, as.character(results$matching[matching]))
    blood_clonotypes <- c(blood_clonotypes, as.character(results$normal_clonotypes[matching]))
  }  
}

