
exome_ref <- read.csv("/Users/Eric/tcrseq/new/reference/cghub_exome_primary_tumorsamples_3-18-16_sorted.tsv", sep="\t")

gsva_dir <- "/Users/Eric/tcrseq/new/gsva/"
clinical_dir <- "/Users/Eric/tcrseq/new/clinical_processed/"
imseq_dir <- "/Users/Eric/tcrseq/old/imseq/"
exome_dir <- "/Users/Eric/tcrseq/old/summaries/"

exome_ref <- exome_ref[!duplicated(exome_ref$participant_id),]
cancer_type_last <- ""
gsva_file <- ""
clinical_file <- ""
results <- c()
for (i in 1:nrow(exome_ref))
{
  cancer_type <- as.character(exome_ref$disease[i])
  if (cancer_type!=cancer_type_last) {
    gsva_path <- paste(gsva_dir, cancer_type, sep="")
    gsva_extension <- paste(gsva_path, "_gsva.txt", sep="")
    gsva_file <- read.csv(gsva_extension, sep="\t")
    clinical_path <- paste(clinical_dir, cancer_type, sep="")
    clinical_extension <- paste(clinical_path, "_clinical.txt", sep="")
    clinical_file <- read.csv(clinical_extension, sep="\t", header=FALSE) 
    rownames(clinical_file) <- clinical_file[,1]
    clinical_file <- clinical_file[,-1]
  }
  
  imseq_reads <- NA
  imseq_clonotypes <- NA
  exome_id <- as.character(exome_ref$analysis_id[i])
  imseq_path <- paste(imseq_dir, exome_id, sep="")
  imseq_extension <- paste(imseq_path, ".tsv", sep="")
  if (file.exists(imseq_extension)) {
    imseq_file <- read.csv(imseq_extension, sep="\t")
    imseq_reads <- nrow(imseq_file)
    if (imseq_reads>0) {
      imseq_clonotypes <- length(unique(imseq_file[,10]))
    } else {
      imseq_clonotypes <- 0
    }
  }
  
  exome_reads <- NA
  exome_path <- paste(exome_dir, exome_id, sep="")
  exome_extension <- paste(exome_path, ".txt", sep="")
  if (file.exists(exome_extension)) {
    exome_file <- read.csv(exome_extension, sep="\t", header=FALSE)
    exome_reads <- sum(c(exome_file[,3], exome_file[,4]))
  }
  
  barcode <- exome_ref$barcode[i]
  barcode <- substr(barcode, 1, 16)
  gsva_barcodes <- gsub(".", "-", colnames(gsva_file), fixed=TRUE)
  gsva_barcodes <- substr(gsva_barcodes,1,16)
  gsva_index <- which(gsva_barcodes==barcode)
  gsva_results <- rep(NA, 22)
  if (length(gsva_index)>0)
  {
    gsva_results <- gsva_file[,gsva_index]
  }
  
  uuid <- as.character(exome_ref$participant_id[i])
  clinical_col <- which(as.matrix(clinical_file[2,])==uuid)
  clinical_results <- rep(NA, 4)
  if (length(clinical_col>0)) {
    clinical_results <- as.numeric(as.matrix(clinical_file[c(3,4,5,6), clinical_col]))
  }
  
  results_row <- c(uuid, barcode, cancer_type, exome_id, exome_reads, imseq_reads, imseq_clonotypes, clinical_results, gsva_results)
  results <- rbind(results, results_row)
}

pancan_colnames <- c("patient_uuid", "sample_barcode", "cohort", "exome_uuid", "exome_reads",
                     "exome_cdr3", "exome_clonotypes", "days_to_birth", "days_to_death",
                     "days_to_last_followup", "percent_tils", "T-regs", "T-CD8d", "T-CD4n", 
                     "T-CD4mr", "T-CD4ma", "T-hf", "T-gd", "B-n", "B-m", "Plasma","NK-r", "NK-a", 
                     "Mono", "Macro-0", "Macro-1", "Macro-2", "Den-r", "Den-a", "Mast-r", "Mast-a",
                      "Eos", "Neut")
colnames(results) <- pancan_colnames

write.table(results, "/Users/Eric/tcrseq/new/processed/pancan_exome_gsva_clinical_dedup.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)

