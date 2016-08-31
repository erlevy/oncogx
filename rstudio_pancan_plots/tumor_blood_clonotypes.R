library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)
library(reshape2)
library(dplyr)
library(ggplot2)
library(plotrix)
library(VennDiagram)
library(ggrepel)

# run clonotypes_analysis.R first

results <- read.csv("/Users/Eric/tcrseq/new/final/pancan_results_8-10-2016.txt", sep="\t")
clonotypes <- read.csv("/Users/Eric/tcrseq/new/final/pancan_clonotypes_8-10-2016.txt", sep="\t")

cdr3_results_current <- results
cancer_current <- cdr3_results_current
# remove cohort size <400
# BLCA, BRCA, COAD, HNSC, KIRC, LGG, LUAD, LUSC, OV, PRAD, STAD, THCA, UCEC
#cancer_current <- filter(cdr3_results_current, cohort=="BLCA" | cohort=="BRCA" | cohort=="COAD" | cohort=="HNSC" |
#                       cohort=="KIRC" | cohort=="LGG" | cohort=="LUAD" | cohort=="LUSC" | cohort=="OV" |
#                       cohort=="PRAD" | cohort=="STAD" | cohort=="THCA" | cohort=="UCEC")
#cancer_current$cohort <- factor(cancer_current$cohort)
#cancer_split <- split(cancer_current, cancer_current$cohort)

# remove THYM since it only has ONE exome case, UVM and CHOL have no RNA
cancer_current <- filter(cdr3_results_current, cohort != "THYM", cohort != "UVM", cohort != "CHOL")
cancer_current$cohort <- factor(cancer_current$cohort)
cancer_split <- split(cancer_current, cancer_current$cohort)

# exome vs. blood
# run clonotypes_analysis.R first
exome_pos <- filter(clonotypes, exome > 0)
blood_pos <- filter(clonotypes, blood > 0)

exome_only <- filter(clonotypes, exome > 0)
exome_pos_n <- length(unique(exome_pos$aa))
blood_pos_n <- length(unique(blood_pos$aa))
both_n <- length(unique(intersect(exome_pos$aa, blood_pos$aa)))

# distribution of cohorts with "sticky" clonotype
shared_patients <- exome_blood_shared_patients_only_c66
shared_patients_results <- filter(cancer_current, patient_uuid %in% shared_patients)
exome_blood_or_n_patients <- filter(cancer_current, !is.na(exome_rpm) | !is.na(blood_rpm))$patient_uuid
exome_blood_baseline <- filter(cancer_current, patient_uuid %in% exome_blood_or_n_patients)
exome_blood_baseline <- mutate(exome_blood_baseline,has_c66=ifelse(patient_uuid %in% shared_patients,"1","0"))
shared_patients_df <- data.frame(table(shared_patients_results$cohort), table(exome_blood_baseline$cohort))
shared_patients_df <- filter(shared_patients_df, Freq>0)
colnames(shared_patients_df)[1] <- "Cohort"
shared_patients_df$Cohort <- factor(shared_patients_df$Cohort)
shared_patients_df$Freq.2 <- shared_patients_df$Freq/shared_patients_df$Freq.1

### FIGURE: c66 cohort fraction barplot ###

ggplot(shared_patients_df, aes(x=reorder(Cohort, -Freq.2), y=Freq.2)) + geom_bar(stat="identity") + labs(x="Cohort", y="Fraction of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

### FIGURE: c66 cohort fraction barplot by gender
shared_patients_gender_df <- data.frame(table(shared_patients_results$cohort, shared_patients_results$gender))
shared_patients_gender_df <- filter(shared_patients_gender_df, Freq>0)
ggplot(shared_patients_gender_df, aes(x=reorder(Var1, -Freq), y=Freq, fill=Var2)) + geom_bar(stat="identity") + labs(x="Cohort", y="Number of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

# distribution of cohorts without "sticky" clonotype
non_shared_patients <- exome_blood_shared_patients_no_c66
non_shared_patients_results <- filter(cancer_current, patient_uuid %in% non_shared_patients)
non_shared_patients_df <- data.frame(table(non_shared_patients_results$cohort), table(exome_blood_baseline$cohort))
non_shared_patients_df <- filter(non_shared_patients_df, Freq>0)
colnames(non_shared_patients_df)[1] <- "Cohort"
non_shared_patients_df$Cohort <- factor(non_shared_patients_df$Cohort)
non_shared_patients_df$Freq.2 <- non_shared_patients_df$Freq/non_shared_patients_df$Freq.1

### FIGURE: non c66 cohort fraction barplot ###
ggplot(non_shared_patients_df, aes(x=reorder(Cohort, -Freq.2), y=Freq.2)) + geom_bar(stat="identity") + labs(x="Cohort", y="Fraction of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

### FIGURE: non c66 cohort fraction barplot by gender
non_shared_patients_gender_df <- data.frame(table(non_shared_patients_results$cohort, non_shared_patients_results$gender))
non_shared_patients_gender_df <- filter(non_shared_patients_gender_df, Freq>0)
ggplot(non_shared_patients_gender_df, aes(x=reorder(Var1, -Freq), y=Freq, fill=Var2)) + geom_bar(stat="identity") + labs(x="Cohort", y="Number of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

# profile gender with c66
shared_patients_gender_df <- data.frame(table(shared_patients_results$gender), table(exome_blood_baseline$gender))
shared_patients_gender_table <- matrix(c(shared_patients_gender_df[1,2], shared_patients_gender_df[1,4],
                                         shared_patients_gender_df[2,2], shared_patients_gender_df[2,4]),
                                       nrow=2, byrow=TRUE)
non_shared_patients_gender_df <- data.frame(table(non_shared_patients_results$gender), table(exome_blood_baseline$gender))
non_shared_patients_gender_table <- matrix(c(non_shared_patients_gender_df[1,2], non_shared_patients_gender_df[1,4],
                                             non_shared_patients_gender_df[2,2], non_shared_patients_gender_df[2,4]), nrow=2, byrow=TRUE)

public_clonotypes_results <- filter(cancer_current, as.character(patient_uuid) %in% public_patients)
exome_rna_blood_or_n_patients <- filter(cancer_current, !is.na(rna_rpm) | !is.na(exome_rpm) | !is.na(blood_rpm))$patient_uuid
public_clonotypes_baseline <- filter(cancer_current, patient_uuid %in% exome_rna_blood_or_n_patients)
table(public_clonotypes_results$cohort)

# coverage vs c66
exome_blood_non_c66 <- filter(cancer_current, patient_uuid %in% exome_or_blood_patients)
exome_blood_non_c66 <- filter(exome_blood_non_c66, !(patient_uuid %in% exome_blood_shared_patients_only_c66))
t.test(shared_patients_results$exome_reads, exome_blood_non_c66$exome_reads)
t.test(shared_patients_results$blood_reads, exome_blood_non_c66$blood_reads)

exome_blood_non_public <- filter(cancer_current, patient_uuid %in% exome_or_blood_patients)
exome_blood_non_public <- filter(exome_blood_non_public, !(patient_uuid %in% exome_blood_shared_patients_no_c66))
t.test(shared_patients_results$exome_reads, exome_blood_non_c66$exome_reads)
t.test(shared_patients_results$blood_reads, exome_blood_non_c66$blood_reads)

### new stuff

# correlation between diversity in blood and tumor exome?
tumor_blood_exome_pos <- filter(results, blood_clonotypes > 0 & exome_clonotypes > 0)
# only cases with >0 clonotypes in both tumor and blood
blood_clonotype_diversity <- tumor_blood_exome_pos$blood_clonotypes/tumor_blood_exome_pos$blood_cdr3
tumor_exome_clonotype_diversity <- tumor_blood_exome_pos$exome_clonotypes/tumor_blood_exome_pos$exome_cdr3


