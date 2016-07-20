
clinical_dir <- "/Users/Eric/tcrseq/new/clinical_raw"

files_clinical <- list.files(path=clinical_dir, pattern="*.txt", full.names=T, recursive=FALSE)

clinical1 <- read.csv(files_clinical[1], sep="\t", header=FALSE)

run_clinical <- function(clinical_path)
{
  clinical1 <- read.csv(clinical_path, sep="\t", header=FALSE)
  # what we need: patient barcode, patient uuid, days to birth, days to death, days to last followup
  clinical_fields <- c("patient.bcr_patient_barcode", "patient.bcr_patient_uuid", "patient.days_to_birth",
                       "patient.days_to_death", "patient.days_to_last_followup", "patient.samples.sample.portions.portion.slides.slide.percent_lymphocyte_infiltration")
  
  clinical_processed <- c()
  for (i in 1:length(clinical_fields))
  {
    field_match <- which(clinical1$V1==clinical_fields[i])
    clinical_processed <- rbind(clinical_processed, clinical1[field_match,])
  }
  
  cancer_type <- basename(strsplit(clinical_path, '[.]')[[1]][1])
  output_dir <- "/Users/Eric/tcrseq/new/clinical_processed/"
  output_file <- paste(cancer_type, "_clinical.txt", sep="")
  clinical_out <- paste(output_dir, output_file, sep="")
  write.table(clinical_processed, clinical_out, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#  return(clinical_processed)
}

lapply(files_clinical, run_clinical)
