
exome_ref <- read.csv("/Users/Eric/tcrseq/new/reference/exome_ref_unique.txt", sep="\t")
rna_ref <- read.csv("/Users/Eric/tcrseq/new/reference/rna_ref_unique.txt", sep="\t")
blood_ref <- read.csv("/Users/Eric/tcrseq/new/reference/blood_ref_unique.txt", sep="\t")

gsva_dir <- "/Users/Eric/tcrseq/new/gsva/"
clinical_dir <- "/Users/Eric/tcrseq/new/clinical_processed/"
exome_imseq_dir <- "/Users/Eric/tcrseq/new/imseq/exome/"
exome_summary_dir <- "/Users/Eric/tcrseq/new/summaries/exome/"
rna_imseq_dir <- "/Users/Eric/tcrseq/new/imseq/rna/"
rna_summary_dir <- "/Users/Eric/tcrseq/new/summaries/rna/"
rna_tcrb_dir <- "/Users/Eric/tcrseq/new/summaries/tcrb/"
blood_imseq_dir <- "/Users/Eric/tcrseq/new/imseq/blood/"
blood_summary_dir <- "/Users/Eric/tcrseq/new/summaries/blood/"

cancer_type_last <- ""
gsva_file <- ""
clinical_file <- ""
results <- matrix(NA, nrow=nrow(rna_ref), ncol=42)
for (i in 1:nrow(rna_ref))
{
  cancer_type <- as.character(rna_ref$disease[i])
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
  
  rna_id <- as.character(rna_ref$analysis_id[i])
  rna_imseq_path <- paste(rna_imseq_dir, rna_id, sep="")
  rna_imseq_extension <- paste(rna_imseq_path, ".tsv", sep="")
  rna_imseq_reads <- NA
  rna_imseq_clonotypes <- NA
  rna_summary_reads <- NA
  rna_tcrb_reads <- NA
  if (file.exists(rna_imseq_extension)) {
    if (file.info(rna_imseq_extension)$size>0)
    {
      rna_imseq_file <- read.csv(rna_imseq_extension, sep="\t")
      rna_imseq_reads <- nrow(rna_imseq_file)
      if (rna_imseq_reads>0) {
        rna_imseq_clonotypes <- length(unique(rna_imseq_file[,10]))
      } else {
        rna_imseq_clonotypes <- 0
      }        
    }
  } 
  rna_summary_path <- paste(rna_summary_dir, rna_id, sep="")
  rna_summary_extension <- paste(rna_summary_path, ".txt", sep="")
  if (file.exists(rna_summary_extension)) {
    if (file.info(rna_summary_extension)$size>0)
    {
      rna_summary_file <- read.csv(rna_summary_extension, sep="\t", header=FALSE)
      rna_summary_reads <- sum(c(rna_summary_file[,3], rna_summary_file[,4]))        
    }
  }
  rna_tcrb_path <- paste(rna_tcrb_dir, rna_id, sep="")
  rna_tcrb_extension <- paste(rna_tcrb_path, ".idxstats.txt", sep="")
  if (file.exists(rna_tcrb_extension)) {
    if (file.info(rna_tcrb_extension)$size>0)
    {
      rna_tcrb_file <- read.csv(rna_tcrb_extension, sep="\t", header=FALSE, row.names=1)
      rna_tcrb_reads <- rna_tcrb_file["chr7",2]       
    }
  }
  
  exome_id <- NA
  exome_imseq_reads <- NA
  exome_imseq_clonotypes <- NA
  exome_summary_reads <- NA
  exome_match <- which(as.character(exome_ref$participant_id)==as.character(rna_ref$participant_id[i]))
  if (length(exome_match)!=0)
  {
    exome_id <- as.character(exome_ref$analysis_id[exome_match])
    exome_imseq_path <- paste(exome_imseq_dir, exome_id, sep="")
    exome_imseq_extension <- paste(exome_imseq_path, ".tsv", sep="")
    if (file.exists(exome_imseq_extension)) {
      if (file.info(exome_imseq_extension)$size>0)
      {
        exome_imseq_file <- read.csv(exome_imseq_extension, sep="\t")
        exome_imseq_reads <- nrow(exome_imseq_file)
        if (exome_imseq_reads>0) {
          exome_imseq_clonotypes <- length(unique(exome_imseq_file[,10]))
        } else {
          exome_imseq_clonotypes <- 0
        }      
      }
    }
    exome_summary_path <- paste(exome_summary_dir, exome_id, sep="")
    exome_summary_extension <- paste(exome_summary_path, ".txt", sep="")
    if (file.exists(exome_summary_extension)) {
      if (file.info(exome_summary_extension)$size>0)
      {
        exome_summary_file <- read.csv(exome_summary_extension, sep="\t", header=FALSE)
        exome_summary_reads <- sum(c(exome_summary_file[,3], exome_summary_file[,4]))      
      }
    }
  }
  
  blood_id <- NA
  blood_imseq_reads <- NA
  blood_imseq_clonotypes <- NA
  blood_summary_reads <- NA
  blood_match <- which(as.character(blood_ref$participant_id)==as.character(rna_ref$participant_id[i]))
  if (length(blood_match)!=0)
  {
    blood_id <- as.character(blood_ref$analysis_id[blood_match])
    blood_imseq_path <- paste(blood_imseq_dir, blood_id, sep="")
    blood_imseq_extension <- paste(blood_imseq_path, ".tsv", sep="")
    if (file.exists(blood_imseq_extension)) {
      if (file.info(blood_imseq_extension)$size>0)
      {
        blood_imseq_file <- read.csv(blood_imseq_extension, sep="\t")
        blood_imseq_reads <- nrow(blood_imseq_file)
        if (blood_imseq_reads>0) {
          blood_imseq_clonotypes <- length(unique(blood_imseq_file[,10]))
        } else {
          blood_imseq_clonotypes <- 0
        }        
      }
    } 
    blood_summary_path <- paste(blood_summary_dir, blood_id, sep="")
    blood_summary_extension <- paste(blood_summary_path, ".txt", sep="")
    if (file.exists(blood_summary_extension)) {
      if (file.info(blood_summary_extension)$size>0)
      {
        blood_summary_file <- read.csv(blood_summary_extension, sep="\t", header=FALSE)
        blood_summary_reads <- sum(c(blood_summary_file[,3], blood_summary_file[,4]))        
      }
    }
  }

  barcode <- rna_ref$barcode[i]
  barcode <- substr(barcode, 1, 16)
  gsva_barcodes <- gsub(".", "-", colnames(gsva_file), fixed=TRUE)
  gsva_barcodes <- substr(gsva_barcodes,1,16)
  gsva_index <- which(gsva_barcodes==barcode)
  gsva_results <- rep(NA, 22)
  if (length(gsva_index)>0)
  {
    gsva_results <- gsva_file[,gsva_index]
  }
  
  uuid <- as.character(rna_ref$participant_id[i])
  clinical_col <- which(as.matrix(clinical_file[2,])==uuid)
  clinical_results <- rep(NA, 4)
  if (length(clinical_col>0)) {
    clinical_results <- as.numeric(as.matrix(clinical_file[c(3,4,5,6), clinical_col]))
  }
  
  results_row <- c(uuid, barcode, cancer_type, exome_id, exome_summary_reads, exome_imseq_reads,
                   exome_imseq_clonotypes, rna_id, rna_summary_reads, rna_imseq_reads, rna_imseq_clonotypes,
                   rna_tcrb_reads, blood_id, blood_summary_reads, blood_imseq_reads, blood_imseq_clonotypes, 
                   clinical_results, gsva_results)
  results[i,] <- results_row
}

pancan_colnames <- c("patient_uuid", "sample_barcode", "cohort", "exome_uuid", "exome_reads",
                     "exome_cdr3", "exome_clonotypes", "rna_uuid", "rna_reads", "rna_cdr3", "rna_clonotypes",
                     "tcrb_reads", "blood_uuid", "blood_reads", "blood_cdr3", "blood_clonotypes",
                     "days_to_birth", "days_to_death","days_to_last_followup", "percent_tils", 
                     "T-regs", "T-CD8d", "T-CD4n", "T-CD4mr", "T-CD4ma", "T-hf", "T-gd", "B-n", 
                     "B-m", "Plasma","NK-r", "NK-a", "Mono", "Macro-0", "Macro-1", "Macro-2", 
                     "Den-r", "Den-a", "Mast-r", "Mast-a","Eos", "Neut")
colnames(results) <- pancan_colnames

write.table(results, "/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-5-2016_RNA.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)

