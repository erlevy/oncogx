
exome_ref <- read.csv("/Users/Eric/tcrseq/new/reference/cghub_exome_primary_tumorsamples_3-18-16_sorted.tsv", sep="\t")
rna_ref <- read.csv("/Users/Eric/tcrseq/new/reference/cghub_rna_primary_tumor_samples_4-5-16.tsv", sep="\t")
blood_ref <- read.csv("/Users/Eric/tcrseq/new/reference/cghub_blood_normal_samples_4-5-16.tsv", sep="\t")

# get 1 exome per sample barcode
unique_patients <- unique(as.character(exome_ref$participant_id))
exome_filtered <- c()
for (i in 1:length(unique_patients))
{
  patient_id <- unique_patients[i]
  patient_rownums <- which(as.character(exome_ref$participant_id)==patient_id)
  patient_rows <- exome_ref[patient_rownums,]
  if (nrow(patient_rows)==1)
  {
    exome_filtered <- rbind(exome_filtered, patient_rows)
  }
  else
  {
    dates <- as.Date(as.character(patient_rows$uploaded), "%m/%d/%y")
    max_dates <- which(dates==max(dates))
    if (length(max_dates)>1)
    {
      exome_filtered <- rbind(exome_filtered, patient_rows[max_dates[1],])         
    }
    else
    {
      exome_filtered <- rbind(exome_filtered, patient_rows[max_dates,])                   
    }
  }
}
write.table(exome_filtered, "exome_ref_unique.txt", quote=FALSE, row.names=FALSE, sep="\t")

# get 1 rna per sample barcode
unique_patients <- unique(as.character(rna_ref$participant_id))
rna_filtered <- c()
for (i in 1:length(unique_patients))
{
  patient_id <- unique_patients[i]
  patient_rownums <- which(as.character(rna_ref$participant_id)==patient_id)
  patient_rows <- rna_ref[patient_rownums,]
  if (nrow(patient_rows)==1)
  {
    rna_filtered <- rbind(rna_filtered, patient_rows)
  }
  else
  {
    dates <- as.Date(as.character(patient_rows$uploaded))
    max_dates <- which(dates==max(dates))
    if (length(max_dates)>1)
    {
      rna_filtered <- rbind(rna_filtered, patient_rows[max_dates[1],])         
    }
    else
    {
      rna_filtered <- rbind(rna_filtered, patient_rows[max_dates,])                   
    }
  }
}
write.table(rna_filtered, "rna_ref_unique.txt", quote=FALSE, row.names=FALSE, sep="\t")


# get 1 blood per sample barcode
unique_patients <- unique(as.character(blood_ref$participant_id))
blood_filtered <- c()
for (i in 1:length(unique_patients))
{
  patient_id <- unique_patients[i]
  patient_rownums <- which(as.character(blood_ref$participant_id)==patient_id)
  patient_rows <- blood_ref[patient_rownums,]
  if (nrow(patient_rows)==1)
  {
    blood_filtered <- rbind(blood_filtered, patient_rows)
  }
  else
  {
    dates <- as.Date(as.character(patient_rows$uploaded))
    max_dates <- which(dates==max(dates))
    if (length(max_dates)>1)
    {
      blood_filtered <- rbind(blood_filtered, patient_rows[max_dates[1],])         
    }
    else
    {
      blood_filtered <- rbind(blood_filtered, patient_rows[max_dates,])                   
    }
  }
}
write.table(blood_filtered, "blood_ref_unique.txt", quote=FALSE, row.names=FALSE, sep="\t")
