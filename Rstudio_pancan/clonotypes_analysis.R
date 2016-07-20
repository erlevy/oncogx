library(dplyr)

clonotypes <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_all_7-14-2016.txt", sep="\t")

clonotypes$clonotype <- as.character(clonotypes$clonotype)
clonotypes$patient_id <- as.character(clonotypes$patient_id)
clonotypes$aa <- as.character(clonotypes$aa)

# compare exome and blood
exome_pos_blood_pos <- filter(clonotypes, blood > 0 & exome > 0)
exome_neg_blood_pos <- filter(clonotypes, blood > 0 & exome == 0)
exome_pos_blood_neg <- filter(clonotypes, blood == 0 & exome > 0)

clonotype_common <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
shared <- clonotypes[which(clonotypes$clonotype==clonotype_common),]

# compare exome and rna
exome_pos_rna_pos <- filter(clonotypes, rna > 0 & exome > 0)
exome_neg_rna_pos <- filter(clonotypes, rna > 0 & exome == 0)
exome_pos_rna_neg <- filter(clonotypes, rna == 0 & exome > 0)

