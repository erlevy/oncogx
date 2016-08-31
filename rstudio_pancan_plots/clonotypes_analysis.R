library(dplyr)
library(ggplot2)
library(colorspace)
library(viridis)

results <- read.csv("/Users/Eric/tcrseq/new/final/pancan_results_8-10-2016.txt", sep="\t")
clonotypes <- read.csv("/Users/Eric/tcrseq/new/final/pancan_clonotypes_8-10-2016.txt", sep="\t")

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

### new paper work start here 8/21 ###

# just tumor exome
exome_pos <- filter(clonotypes, exome > 0)
tumor_exome_cl <- unique(exome_pos$clonotype)
length(tumor_exome_cl)
length(unique(exome_pos$patient_id))

tumor_exome_shared_aa <- names(which(table(exome_pos$clonotype)>1))
exome_tumor_shared <- filter(exome_pos, clonotype %in% tumor_exome_shared_aa)
length(unique(exome_tumor_shared$clonotype))
length(unique(exome_tumor_shared$patient_id))
sort(table(exome_tumor_shared$clonotype))

hist((sort(table(exome_tumor_shared$clonotype), decreasing = TRUE))) # 5 AA one is likely artifact
ggplot(melt(table(exome_tumor_shared$clonotype)), aes(value)) + geom_histogram() + labs(x="Number of patients with clonotype", y="Number of clonotypes") + theme(text=element_text(size=16))

# just tumor RNA
rna_pos <- filter(clonotypes, rna > 0)
tumor_rna_cl <- unique(rna_pos$clonotype)
length(tumor_rna_cl)
length(unique(rna_pos$patient_id))

tumor_rna_shared_aa <- names(which(table(rna_pos$clonotype)>1))
rna_tumor_shared <- filter(rna_pos, clonotype %in% tumor_rna_shared_aa)
length(unique(rna_tumor_shared$clonotype))
length(unique(rna_tumor_shared$patient_id))
sort(table(rna_tumor_shared$clonotype))

hist((sort(table(rna_tumor_shared$clonotype), decreasing = TRUE)[-1])) # 5 AA one is likely artifact
ggplot(melt(table(rna_tumor_shared$clonotype)), aes(value)) + geom_histogram() + labs(x="Number of patients with clonotype", y="Number of clonotypes") + theme(text=element_text(size=16))

# tumor exome and RNA
tumor_exome_rna_shared <- intersect(tumor_exome_cl, tumor_rna_cl)
length(tumor_exome_rna_shared)

exome_rna_tumor_shared <- filter(clonotypes, clonotype %in% tumor_exome_rna_shared)
length(unique(exome_rna_tumor_shared$clonotype))
length(unique(exome_rna_tumor_shared$patient_id))
sort(table(exome_rna_tumor_shared$clonotype))

hist((sort(table(exome_rna_tumor_shared$clonotype), decreasing = TRUE)[-1])) # 5 AA one is likely artifact
ggplot(melt(table(exome_rna_tumor_shared$clonotype)), aes(value)) + geom_histogram() + labs(x="Number of patients with clonotype", y="Number of clonotypes") + theme(text=element_text(size=16))

# profile these shared cases
exome_rna_tumor_shared_patients <- unique(exome_rna_tumor_shared$patient_id)
exome_rna_tumor_shared_data <- filter(results, patient_uuid %in% exome_rna_tumor_shared_patients)

overall_cohorts <- table(results$cohort)
shared_cohorts <- table(exome_rna_tumor_shared_data$cohort)
shared_df <- data.frame(overall_cohorts-shared_cohorts, shared_cohorts)
shared_df <- filter(shared_df, Freq.1>0)
shared_df$Var1.1 <- NULL
colnames(shared_df) <- c("Var1", "Non-shared", "Shared")
shared_df_melt <- melt(shared_df)

ggplot(shared_df_melt, aes(x=Var1, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + labs(x="Cohort", y="Fraction") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

# just blood exome
blood_pos <- filter(clonotypes, blood > 0)
blood_exome_cl <- unique(blood_pos$clonotype)
length(blood_exome_cl)
length(unique(blood_pos$patient_id))

blood_exome_shared_aa <- names(which(table(blood_pos$clonotype)>1))
exome_blood_shared <- filter(blood_pos, clonotype %in% blood_exome_shared_aa)
length(unique(exome_blood_shared$clonotype))
length(unique(exome_blood_shared$patient_id))
sort(table(exome_blood_shared$clonotype))

hist((sort(table(exome_blood_shared$clonotype), decreasing = TRUE))) # 5 AA one is likely artifact
ggplot(melt(table(exome_blood_shared$clonotype)), aes(value)) + geom_histogram() + labs(x="Number of patients with clonotype", y="Number of clonotypes") + theme(text=element_text(size=16))

# tumor exome and blood
tumor_exome_blood_shared <- intersect(tumor_exome_cl, blood_exome_cl)
length(tumor_exome_blood_shared)

exome_blood_tumor_shared <- filter(clonotypes, clonotype %in% tumor_exome_blood_shared)
length(unique(exome_blood_tumor_shared$clonotype))
length(unique(exome_blood_tumor_shared$patient_id))
sort(table(exome_blood_tumor_shared$clonotype))

hist((sort(table(exome_blood_tumor_shared$clonotype), decreasing = TRUE)[-1])) # 5 AA one is likely artifact
ggplot(melt(table(exome_blood_tumor_shared$clonotype)), aes(value)) + geom_histogram() + labs(x="Number of patients with clonotype", y="Number of clonotypes") + theme(text=element_text(size=16))

# all three
tumor_all_shared <- intersect(tumor_exome_blood_shared, tumor_exome_rna_shared)
length(tumor_all_shared)

all_shared <- filter(clonotypes, clonotype %in% tumor_all_shared)
length(unique(all_shared$clonotype))
length(unique(all_shared$patient_id))
sort(table(all_shared$clonotype))

hist((sort(table(all_shared$clonotype), decreasing = TRUE)[-1])) # 5 AA one is likely artifact
ggplot(melt(table(all_shared$clonotype)), aes(value)) + geom_histogram() + labs(x="Number of patients with clonotype", y="Number of clonotypes") + theme(text=element_text(size=16))

# tumor exome and blood
blood_pos <- filter(clonotypes, blood > 0)
tumor_blood_cl <- unique(blood_pos$clonotype)
length(tumor_blood_cl)
length(unique(blood_pos$patient_id))

tumor_all_shared <- intersect(tumor_exome_cl, tumor_blood_cl)
length(tumor_all_shared)

exome_blood_tumor_shared <- filter(clonotypes, clonotype %in% tumor_all_shared)
length(unique(exome_blood_tumor_shared$clonotype))
length(unique(exome_blood_tumor_shared$patient_id))
sort(table(exome_blood_tumor_shared$clonotype))

# tumor exome and blood c66 comparison
clonotype_common <- "CATSRDTELQCFLLSVHKPHCFPDPGAFS"

# clonal diversity by cohort
exome_pos_results <- filter(results, exome_clonotypes>0)
rna_pos_results <- filter(results, rna_clonotypes>0)
ggplot(exome_pos_results, aes(exome_cdr3, exome_clonotypes, colour=cohort)) + geom_point() + labs(x="exome CDR3 reads", y="Number of clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(rna_cdr3, rna_clonotypes, colour=cohort)) + geom_point() + labs(x="RNA CDR3 reads", y="Number of clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

ggplot(rna_pos_results, aes(log10(tcrb_rpm), log10(rna_rpm), colour=cohort)) + geom_point() + labs(x="TCRB RNA rpm (log10)", y="RNA CDR3 rpm (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_cdr3), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA CDR3 reads (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_clonotypes), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA CDR3 clonotypes (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

cor.test(exome_pos_results$exome_reads, exome_pos_results$exome_clonotypes)
cor.test(rna_pos_results$tcrb_rpm, rna_pos_results$rna_rpm)
cor.test(rna_pos_results$tcrb_reads, rna_pos_results$rna_clonotypes)
