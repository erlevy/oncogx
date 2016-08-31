library(dplyr)

burden_dir <- "/Users/Eric/tcrseq/new/mutation/mutation_sample/"
results <- read.csv("/Users/Eric/tcrseq/new/final/pancan_results_8-10-2016.txt", sep="\t")

results <- arrange(results, cohort)

burden <- rep(NA, nrow(results))
cancer_type_last <- ""
mutation_file <- c()
for (i in 1:nrow(results))
{
  cancer_type <- as.character(results$cohort[i])
  if (cancer_type!=cancer_type_last) {
    mutation_start <- paste(burden_dir, cancer_type, sep="")
    mutation_extension <- paste(mutation_start, "-TP.samplefeatures.txt", sep="")
    if (file.exists(mutation_extension)) {
      mutation_file <- read.csv(mutation_extension, sep="\t")
    }
  }
  
  nonsyn <- NA
  patient_barcode <- substr(as.character(results$sample_barcode[i]), 1, 12)
  # rate_non
  row_match <- which(as.character(mutation_file[,1])==patient_barcode)
  
  if (length(row_match)>0) {nonsyn <- mutation_file[row_match, "rate_non"]}
  burden[i] <- nonsyn
  
  if ((i %% 1000)==0) {print(i)}
}

results_burden <- cbind(results, burden)
colnames(results_burden)[57] <- "rate_non"

write.table(results_burden, "/Users/Eric/tcrseq/new/final/pancan_results_8-28-2016.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)

