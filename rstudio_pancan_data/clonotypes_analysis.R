library(dplyr)

clonotypes <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_all_7-30-2016.txt", sep="\t")

clonotypes$clonotype <- as.character(clonotypes$clonotype)
clonotypes$patient_id <- as.character(clonotypes$patient_id)
clonotypes$nuc <- as.character(clonotypes$nuc)

# compare exome and blood
# defining "public clonotype" as any that is found in both exome and blood (though not in same patient)
#exome_pos_blood_pos <- filter(clonotypes, blood > 0 & exome > 0)
#exome_neg_blood_pos <- filter(clonotypes, blood > 0 & exome == 0)
#exome_pos_blood_neg <- filter(clonotypes, blood == 0 & exome > 0)

#clonotype_common <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
clonotype_common <- "CATSRDTELQCFLLSVHKPHCFPDPGAFS"
shared <- clonotypes[which(clonotypes$clonotype==clonotype_common),]

#public_clonotypes_patients <- unique(exome_pos_blood_pos$patient_id)

exome_pos <- filter(clonotypes, exome > 0)
blood_pos <- filter(clonotypes, blood > 0)

exome_blood_shared_aa <- unique(intersect(exome_pos$clonotype, blood_pos$clonotype))
exome_or_blood <- filter(clonotypes, exome > 0 | blood > 0)
exome_or_blood_patients <- unique(exome_or_blood$patient_id)
exome_blood_public_aa <- names(which(table(exome_or_blood$clonotype)>1))
exome_blood_shared_all <- filter(exome_or_blood, clonotype %in% exome_blood_shared_aa)
exome_blood_shared_patients <- unique(exome_blood_shared_all$patient_id)
exome_blood_public_all <- filter(exome_or_blood, clonotype %in% exome_blood_public_aa)
exome_blood_public_patients <- unique(exome_blood_public_all$patient_id)

# 100 shared aa sequences between blood and tumor exome
# 5934 patients with a clonotype in exome OR blood
# 574 that have a public clonotype
# 348 of these have the "sticky" clonotype

# remove "sticky" clonotype
exome_blood_shared_aa_no_c66 <- exome_blood_shared_aa[exome_blood_shared_aa!=clonotype_common]
exome_blood_shared_no_c66 <- filter(exome_or_blood, clonotype %in% exome_blood_shared_aa_no_c66)
exome_blood_shared_patients_no_c66 <- unique(exome_blood_shared_no_c66$patient_id)

# only "sticky" clonotype
exome_blood_shared_only_c66 <- filter(exome_or_blood, clonotype==clonotype_common)
exome_blood_shared_patients_only_c66 <- unique(exome_blood_shared_only_c66$patient_id)

# includes RNA
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

