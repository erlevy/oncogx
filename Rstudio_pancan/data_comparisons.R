library(stringdist)

### put together both RNA and DNA
results_exome <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-5-2016_exome.txt", sep="\t")
results_RNA <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-5-2016_RNA.txt", sep="\t")

patients_exome_to_add <- setdiff(as.character(results_exome$patient_uuid), as.character(results_RNA$patient_uuid))
rownames(results_exome) <- as.character(results_exome$patient_uuid)
rows_exome_to_add <- results_exome[patients_exome_to_add,]
rownames(rows_exome_to_add) <- NULL
results_both <- rbind(results_RNA, rows_exome_to_add)
write.table(results, "/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-6-2016.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)

### compare
results_new <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-13-2016.txt", sep="\t")

#results_old <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_gsva_clinical_dedup.txt", sep="\t")
results_old <- read.csv("/Users/Eric/BRCA/cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_clustering_four_groups_1078.txt", sep="\t")

results_holt <- read.csv("/Users/Eric/tcrseq/new/holt/TCGA_mitcr_cdr3_result_160401.tsv.datafreeze1.3.filtered.tsv", sep="\t")

clonotypes_new <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_all_7-14-2016.txt", sep="\t")

### comparing to BRCA
comparison <- c()
comparison_ids <- c()
tcr <- c()
for (i in 1:nrow(results_old))
{
  patient_id_old <- as.character(results_old$patient_id[i])
  patient_id_match <- which(as.character(results_new$patient_uuid)==patient_id_old)
  if (length(patient_id_match==1) && !is.na(results_old$exome_imseq[i]) 
      && !is.na(results_old$rna_imseq[i]) && !is.na(results_old$blood_cdr3[i]))
  {
    old_results <- results_old[i,c("exome_reads", "exome_imseq", "exome_clonotypes",
                                   "rna_reads", "rna_imseq", "rna_clonotypes", "tcrb_reads",
                                   "blood_reads", "blood_cdr3", "blood_clonotypes")]
    new_results <- results_new[patient_id_match,c("exome_reads", "exome_cdr3", "exome_clonotypes",
                                                "rna_reads", "rna_cdr3", "rna_clonotypes", "tcrb_reads",
                                                "blood_reads", "blood_cdr3", "blood_clonotypes")]
    subtraction <- as.matrix(old_results-new_results)
    comparison <- rbind(comparison, subtraction)
    comparison_ids <- c(comparison_ids, patient_id_old)
    tcr <- rbind(tcr, c(old_results$tcrb_reads, new_results$tcrb_reads, old_results$rna_imseq))
  }
}

comparison_ids[which(comparison[,2]!=0)]

### comparing to Holt

clonotypes_rna <- clonotypes_new[clonotypes_new$rna>0,]
pancan_rna_unique_nuc <- as.character(unique(clonotypes_rna$clonotype))
holt_unique_nuc <- as.character(unique(results_holt$nucSeq))
length(intersect(pancan_rna_unique_nuc, holt_unique_nuc))

pancan_rna_unique_aa <- as.character(unique(clonotypes_rna$aa))
pancan_rna_unique_aa <- pancan_rna_unique_aa[complete.cases((pancan_rna_unique_aa))]
holt_unique_aa <- as.character(unique(results_holt$aaSeq))
length(intersect(pancan_rna_unique_aa, holt_unique_aa))

# distances_mat <- matrix(NA, nrow=length(pancan_rna_unique_aa), ncol=3)
# for (i in 1:length(pancan_rna_unique_aa))
# {
#   distances <- stringdist(pancan_rna_unique_aa[i], holt_unique_aa)
#   distance <- min(distances)
#   aa_holt <- holt_unique_aa[which(distances==distance)][1]
#   aa_pancan <- pancan_rna_unique_aa[i]
#   distances_mat[i,] <- c(distance, aa_pancan, aa_holt)
#}

#write.table(distances_mat, "pancan_clonotypes_holt_edit_distances.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

distances_mat <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_holt_edit_distances.txt", sep="\t", header=FALSE)

distances <- distances_mat[,1]

hist(distances, xlab="Minimum string distance", main="Distances with Holt")

# within distances
distances_mat <- matrix(NA, nrow=length(pancan_rna_unique_aa), ncol=3)
for (i in 1:length(pancan_rna_unique_aa))
{
 pancan_rna_unique_without <- pancan_rna_unique_aa[-i]
 distances <- stringdist(pancan_rna_unique_aa[i], pancan_rna_unique_without)
 distance <- min(distances)
 aa_without <- pancan_rna_unique_without[which(distances==distance)][1]
 aa_pancan <- pancan_rna_unique_aa[i]
 distances_mat[i,] <- c(distance, aa_pancan, aa_without)
}
 
 
write.table(distances_mat, "pancan_clonotypes_within_edit_distances.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

distances_mat <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_within_edit_distances.txt", sep="\t", header=FALSE)

distances <- distances_mat[,1]

hist(distances, xlab="Minimum string distance", main="Distances within pancan")

# compare to Li
li_file <- "/Users/Eric/tcrseq/new/li/Supplementary_Table_1.fa"
con <- file(li_file, open="r")
lines <- readLines(con)
close(con, type="r")

pancan_rna_aa_c_f <- c()
for (i in 1:length(pancan_rna_unique_aa))
{
  aa <- pancan_rna_unique_aa[i]
  first <- regexpr("C", aa)
  last <- regexpr("F", aa)
  if (first != -1 && last != -1 && first<last) {
    aa_c_f <- substr(aa, first, last)
    pancan_rna_aa_c_f <- c(pancan_rna_aa_c_f, aa_c_f)
  }
}

li_table <- matrix(data=NA, nrow=length(lines)/2, ncol=2)
header <- c()
for (i in 1:length(lines))
{
  line <- lines[i]
  if (substr(line, 1, 1) == ">") {
    header <- line
  } else {
    li_table[i/2,] <- c(header, line)
  }
}

li_aa <- unique(li_table[,2])
li_aa_c_f <- c()
for (i in 1:length(li_aa))
{
  aa <- li_aa[i]
  first <- regexpr("C", aa)
  last <- regexpr("F", aa)
  if (first != -1 && last != -1 && first<last) {
    aa_c_f <- substr(aa, first, last)
    li_aa_c_f <- c(li_aa_c_f, aa_c_f)
    }
}
pancan_rna_aa_c_f_unique <- unique(pancan_rna_aa_c_f)
li_aa_c_f_unique <- unique(li_aa_c_f)
length(intersect(pancan_rna_aa_c_f_unique, li_aa_c_f_unique))

distances_mat <- matrix(NA, nrow=length(pancan_rna_aa_c_f_unique), ncol=3)
for (i in 1:length(pancan_rna_aa_c_f_unique))
{
  distances <- stringdist(pancan_rna_aa_c_f_unique[i], li_aa_c_f_unique)
  distance <- min(distances)
  aa_li <- li_aa_c_f_unique[which(distances==distance)][1]
  aa_pancan <- pancan_rna_aa_c_f_unique[i]
  distances_mat[i,] <- c(distance, aa_pancan, aa_li)
}

write.table(distances_mat, "/Users/Eric/tcrseq/new/processed/pancan_clonotypes_li_edit_distances.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

distances_mat <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_li_edit_distances.txt", sep="\t", header=FALSE)

hist(distances_mat[,1], main="Histogram of string distances", xlab="String distance")
  
  
  