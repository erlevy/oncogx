
ccle_ref <- read.csv("/Users/Eric/tcrseq/new/reference/cghub_ccle_6-29.tsv", sep="\t")

ccle_imseq_dir <- "/Users/Eric/tcrseq/new/ccle/imseq/"
ccle_summary_dir <- "/Users/Eric/tcrseq/new/ccle/summaries/"

results <- matrix(NA, nrow=nrow(ccle_ref), ncol=6)

for (i in 1:nrow(ccle_ref))
{
  ccle_cohort <- as.character(ccle_ref$disease[i])
  ccle_imseq_reads <- NA
  ccle_imseq_clonotypes <- NA
  ccle_id <- as.character(ccle_ref$analysis_id[i])
  ccle_analysis <- as.character(ccle_ref$analyte_type[i])

  ccle_imseq_path <- paste(ccle_imseq_dir, ccle_id, sep="")
  ccle_imseq_extension <- paste(ccle_imseq_path, ".tsv", sep="")
  if (file.exists(ccle_imseq_extension)) {
    if (file.info(ccle_imseq_extension)$size>0)
    {
      ccle_imseq_file <- read.csv(ccle_imseq_extension, sep="\t")
      ccle_imseq_reads <- nrow(ccle_imseq_file)
      if (ccle_imseq_reads>0) {
        ccle_imseq_clonotypes <- length(unique(ccle_imseq_file[,10]))
      } else {
        ccle_imseq_clonotypes <- 0
      }      
    }
  }
  
  ccle_summary_reads <- NA
  ccle_summary_path <- paste(ccle_summary_dir, ccle_id, sep="")
  ccle_summary_extension <- paste(ccle_summary_path, ".txt", sep="")
  if (file.exists(ccle_summary_extension)) {
    if (file.info(ccle_summary_extension)$size>0)
    {
      ccle_summary_file <- read.csv(ccle_summary_extension, sep="\t", header=FALSE)
      ccle_summary_reads <- sum(c(ccle_summary_file[,3], ccle_summary_file[,4]))      
    }
  }
  results_row <- c(ccle_id, ccle_analysis, ccle_cohort, ccle_imseq_reads, ccle_imseq_clonotypes, ccle_summary_reads)
  results[i,] <- results_row
}

ccle_colnames <- c("analysis_id", "analysis_type", "disease", "imseq", "clonotypes", "reads")
colnames(results) <- ccle_colnames

write.table(results, "/Users/Eric/tcrseq/new/processed/ccle_imseq_results_7-22-2016.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)


### get output of all clonotypes by patient id
results <- read.csv("/Users/Eric/tcrseq/new/processed/ccle_imseq_results_7-22-2016.txt", sep="\t")
file_ccle <- list.files(path=ccle_imseq_dir, pattern="*.tsv", full.names=T, recursive=FALSE)
clonotypes_all <- matrix(nrow=0, ncol=4)
colnames(clonotypes_all) <- c("clonotype", "patient_id", "ccle", "aa")
for (i in 1:nrow(results))
{
  clonotypes_patient <- matrix(nrow=0, ncol=4)
  ccle_reads <- matrix(nrow=0, ncol=11)
  ccle_file <- paste(results$analysis_id[i], ".tsv", sep="")
  ccle_path <- paste(ccle_imseq_dir, ccle_file, sep="")
  
  if(!file.exists(ccle_path)) {next}
  ccle_reads <- as.matrix(read.csv(ccle_path, sep="\t", header=T))
  
  if (nrow(ccle_reads)>0) 
  {
    for (j in 1:nrow(ccle_reads))
    {
      cl <- ccle_reads[j,10]
      aa <- ccle_reads[j,11]
      match <- which(clonotypes_patient[,1]==cl)
      left_match <- ccle_reads[j,4]
      right_match <- ccle_reads[j,7]
      if (length(match)==0) 
      {
        clonotypes_patient <- rbind(clonotypes_patient, c(cl, as.character(results$analysis_id[i]), 1, aa))
      } else 
      {
        clonotypes_patient[match,3] <- as.character(as.numeric(clonotypes_patient[match,3]) + 1)
      }
    }    
  }
  clonotypes_all <- rbind(clonotypes_all, clonotypes_patient)
}

write.table(clonotypes_all, "/Users/Eric/tcrseq/new/processed/ccle_clonotypes.txt", sep="\t", quote=FALSE, row.names=FALSE)
