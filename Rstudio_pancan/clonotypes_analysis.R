library(dplyr)

clonotypes <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_all_7-26-2016.txt", sep="\t")

clonotypes$clonotype <- as.character(clonotypes$clonotype)
clonotypes$patient_id <- as.character(clonotypes$patient_id)
clonotypes$nuc <- as.character(clonotypes$nuc)

# compare exome and blood
exome_pos_blood_pos <- filter(clonotypes, blood > 0 & exome > 0)
exome_neg_blood_pos <- filter(clonotypes, blood > 0 & exome == 0)
exome_pos_blood_neg <- filter(clonotypes, blood == 0 & exome > 0)

clonotype_common <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
shared <- clonotypes[which(clonotypes$clonotype==clonotype_common),]

public_clonotypes_patients <- unique(exome_pos_blood_pos$patient_id)

exome_pos <- filter(clonotypes, exome > 0)
blood_pos <- filter(clonotypes, blood > 0)

exome_blood_shared_aa <- unique(intersect(exome_pos$aa, blood_pos$aa))
exome_or_blood <- filter(clonotypes, exome > 0 | blood > 0)
exome_or_blood_patients <- exome_or_blood$patient_id
exome_blood_public_aa <- names(which(table(exome_or_blood$aa)>1))
exome_blood_shared_all <- filter(exome_or_blood, aa %in% exome_blood_shared_aa)
exome_blood_shared_patients <- exome_blood_shared_all$patient_id
exome_blood_public_all <- filter(exome_or_blood, aa %in% exome_blood_public_aa)
exome_blood_public_patients <- unique(exome_blood_public_all$patient_id)

public_aa <- names(which(table(clonotypes$clonotype)>1))
public_aa_n <- c()
for (i in 1:length(public_aa))
{
  public_aa_n <- c(public_aa_n, nrow(filter(clonotypes, clonotype==public_aa[i])))
}

public_all <- filter(clonotypes, clonotype %in% public_aa)
public_patients <- unique(public_all$patient_id)

# compare exome and rna
exome_pos_rna_pos <- filter(clonotypes, rna > 0 & exome > 0)
exome_neg_rna_pos <- filter(clonotypes, rna > 0 & exome == 0)
exome_pos_rna_neg <- filter(clonotypes, rna == 0 & exome > 0)

