library(dplyr)
library(ggplot2)
library(colorspace)
library(viridis)
library(grDevices)

results <- read.csv("/Users/Eric/tcrseq/new/final/pancan_results_8-10-2016.txt", sep="\t")
clonotypes <- read.csv("/Users/Eric/tcrseq/new/final/pancan_clonotypes_8-10-2016.txt", sep="\t")

clonotypes$clonotype <- as.character(clonotypes$clonotype)
clonotypes$patient_id <- as.character(clonotypes$patient_id)
clonotypes$nuc <- as.character(clonotypes$nuc)

cancer_current <- filter(results, cohort != "THYM", cohort != "UVM", cohort != "CHOL",
                         cohort != "ACC", cohort != "CESC", cohort != "DLBC", cohort != "ESCA", 
                         cohort != "KICH", cohort != "LGG", cohort != "LIHC", cohort != "MESO",
                         cohort != "PCPG", cohort != "TGCT", cohort != "UCS")
cancer_current$cohort <- factor(cancer_current$cohort)

results <- cancer_current

# clonal diversity by cohort
exome_pos_results <- filter(results, exome_clonotypes>0)
rna_pos_results <- filter(results, rna_clonotypes>0)
ggplot(exome_pos_results, aes(exome_cdr3, exome_clonotypes, colour=cohort)) + geom_point() + labs(x="Tumor DNA CDR3 reads", y="Number of clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(rna_cdr3, rna_clonotypes, colour=cohort)) + geom_point() + labs(x="Tumor RNA CDR3 reads", y="Number of clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

ggplot(rna_pos_results, aes(log10(tcrb_rpm), log10(rna_rpm), colour=cohort)) + geom_point() + labs(x="TCRB RNA rpm (log10)", y="RNA CDR3 rpm (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_cdr3), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA CDR3 reads (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_clonotypes), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA CDR3 clonotypes (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_clonotypes/rna_rpm), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA clonal diversity (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

exome_rna_pos_results <- filter(results, exome_clonotypes>0, rna_clonotypes>0)
ggplot(exome_rna_pos_results, aes(log10(exome_clonotypes/exome_rpm), log10(rna_clonotypes/rna_rpm), colour=cohort)) + geom_point() + labs(x="Exome clonal diversity (log10)", y="RNA clonal diversity (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

grays <- gray.colors(length(unique(exome_rna_pos_results$cohort)))
grays[16] <- "yellow" # should be STAD
grays[18] <- "red" # should be UCEC
grays[11] <- "blue" # should be PAAD
ggplot(exome_rna_pos_results, aes(log10(exome_clonotypes/exome_rpm), log10(rna_clonotypes/rna_rpm), colour=cohort)) + geom_point() + labs(x="Exome clonal diversity (log10)", y="RNA clonal diversity (log10)") + theme(text=element_text(size=16)) + scale_color_manual(values=grays)

shared_id <- unique(filter(clonotypes, rna > 0 & exome > 0)$patient_id)
shared_id_results <- filter(exome_rna_pos_results, patient_uuid %in% shared_id)
non_shared_id_results <- filter(exome_rna_pos_results, !(patient_uuid %in% shared_id))
shared_id_results$shared <- "shared"
non_shared_id_results$shared <- "no shared"
shared_annotation_results <- rbind(non_shared_id_results, shared_id_results)
shared_annotation_results$shared <- factor(shared_annotation_results$shared, levels = c("shared","no shared"))
shared_annotation_colors <- c("red", "lightgray")
ggplot(shared_annotation_results, aes(log10(exome_clonotypes/exome_rpm), log10(rna_clonotypes/rna_rpm), colour=shared)) + geom_point(alpha = 0.9) + labs(x="Exome clonal diversity (log10)", y="RNA clonal diversity (log10)") + theme(text=element_text(size=16)) + scale_color_manual(values=shared_annotation_colors)

ggplot(exome_rna_pos_results, aes(log10(exome_reads), exome_clonotypes, colour=cohort)) + geom_point() + labs(x="Tumor DNA read count (log10)", y="Tumor DNA clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(exome_rna_pos_results, aes(log10(rna_reads), rna_clonotypes, colour=cohort)) + geom_point() + labs(x="Tumor RNA read count (log10)", y="Tumor RNA clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

cor.test(exome_pos_results$exome_reads, exome_pos_results$exome_clonotypes)
cor.test(rna_pos_results$tcrb_rpm, rna_pos_results$rna_rpm)
cor.test(rna_pos_results$tcrb_reads, rna_pos_results$rna_clonotypes)

# linear regression for # clonotypes by # reads
mylm <- lm(exome_clonotypes~exome_cdr3, data=exome_pos_results)
mylm <- lm(rna_clonotypes~rna_cdr3, data=rna_pos_results)

# dna rna shared clonality statistic
length(shared_id) # 196 with shared DNA RNA clonotype
nrow(exome_pos_results) # 2906 with DNA
nrow(rna_pos_results) # 5699 with RNA
nrow(filter(cancer_current, exome_cdr3 > 0 & rna_cdr3 > 0)) # 2502 with both
length(unique(as.character(filter(clonotypes, exome > 0)$clonotype))) # 5434 exome clonotypes
length(unique(as.character(filter(clonotypes, rna > 0)$clonotype))) # 61914 rna clonotypes

exome_rna_both_results <- filter(cancer_current, exome_cdr3 > 0 & rna_cdr3 > 0)
exome_rna_both_patients <- as.character(exome_rna_both_results$patient_uuid)
exome_rna_both_clonotypes <- filter(clonotypes, patient_id %in% exome_rna_both_patients)
length(unique(as.character(exome_rna_both_clonotypes$clonotype))) # 45996 clonotypes in patients with both exome and rna
length(unique(as.character(filter(exome_rna_both_clonotypes, exome > 0)$clonotype))) # 4616 exome clonotypes in patients with both exome and rna
length(unique(as.character(filter(exome_rna_both_clonotypes, rna > 0)$clonotype))) # 37859 exome clonotypes in patients with both exome and rna

summary(exome_rna_both_results$exome_clonotypes) # mean of 2
summary(exome_rna_both_results$rna_clonotypes) # mean of 16

# check to make sure shared ones don't have different clonal diversity
wilcox.test(shared_id_results$exome_clonotypes/shared_id_results$exome_rpm, non_shared_id_results$exome_clonotypes/non_shared_id_results$exome_rpm)
wilcox.test(shared_id_results$rna_clonotypes/shared_id_results$rna_rpm, non_shared_id_results$rna_clonotypes/non_shared_id_results$rna_rpm)



