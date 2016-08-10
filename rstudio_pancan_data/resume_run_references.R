
ref_complete <- read.csv("/Users/Eric/tcrseq/new/reference/pancancer_rna_ref_2.txt", sep="\t", header=FALSE)

ls_partial <- read.csv("/Users/Eric/tcrseq/new/reference/second/imseq_done_rna_2_processed.txt", sep="\t", header=FALSE)

exome_ref_new <- c()
exomes_complete <- as.character(ls_partial[,9])
exomes_complete <- substr(exomes_complete, 1, 36)

for (i in 1:nrow(ref_complete))
{
  exome_id <- as.character(ref_complete$V1[i])
  if (length(which(exome_id==exomes_complete))==0) {exome_ref_new <- rbind(exome_ref_new, ref_complete[i,])}
}

#write.table(exome_ref_new, "imseq_todo_rna_2.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# starting from the processed data
results_new <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_6-27-2016.txt", sep="\t")
exome_ref <- read.csv("/Users/Eric/tcrseq/new/reference/exome_ref_unique.txt", sep="\t")
rna_ref <- read.csv("/Users/Eric/tcrseq/new/reference/rna_ref_unique.txt", sep="\t")
blood_ref <- read.csv("/Users/Eric/tcrseq/new/reference/blood_ref_unique.txt", sep="\t")

exome_done <- results_new[which(!is.na(results_new$exome_cdr3)),]
rna_done <- results_new[which(!is.na(results_new$rna_cdr3)),]
blood_done <- results_new[which(!is.na(results_new$blood_cdr3)),]

exome_done_uuid <- as.character(exome_done$exome_uuid)
rna_done_uuid <- as.character(rna_done$rna_uuid)
blood_done_uuid <- as.character(blood_done$blood_uuid)

exome_ref_uuid <- as.character(exome_ref$analysis_id)
rna_ref_uuid <- as.character(rna_ref$analysis_id)
blood_ref_uuid <- as.character(blood_ref$analysis_id)

exome_uuid_todo <- setdiff(exome_ref_uuid, exome_done_uuid)
rna_uuid_todo <- setdiff(rna_ref_uuid, rna_done_uuid)
blood_uuid_todo <- setdiff(blood_ref_uuid, blood_done_uuid)

exome_ref_out <- c()
for (i in 1:length(exome_uuid_todo))
{
  ref <- as.character(exome_ref$assembly[which(exome_ref$analysis_id==exome_uuid_todo[i])])
  exome_ref_out <- rbind(exome_ref_out, c(exome_uuid_todo[i], ref))
}

rna_ref_out <- c()
for (i in 1:length(rna_uuid_todo))
{
  ref <- as.character(rna_ref$assembly[which(rna_ref$analysis_id==rna_uuid_todo[i])])
  rna_ref_out <- rbind(rna_ref_out, c(rna_uuid_todo[i], ref))
}

blood_ref_out <- c()
for (i in 1:length(blood_uuid_todo))
{
  ref <- as.character(blood_ref$assembly[which(blood_ref$analysis_id==blood_uuid_todo[i])])
  blood_ref_out <- rbind(blood_ref_out, c(blood_uuid_todo[i], ref))
}

write.table(exome_ref_out, "exome_uuid_todo_final.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(rna_ref_out, "rna_uuid_todo_final.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(blood_ref_out, "blood_uuid_todo_final.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

