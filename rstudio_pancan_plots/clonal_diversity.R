library(dplyr)
library(ggplot2)
library(colorspace)
library(viridis)

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
ggplot(exome_pos_results, aes(exome_cdr3, exome_clonotypes, colour=cohort)) + geom_point() + labs(x="exome CDR3 reads", y="Number of clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(rna_cdr3, rna_clonotypes, colour=cohort)) + geom_point() + labs(x="RNA CDR3 reads", y="Number of clonotypes") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

ggplot(rna_pos_results, aes(log10(tcrb_rpm), log10(rna_rpm), colour=cohort)) + geom_point() + labs(x="TCRB RNA rpm (log10)", y="RNA CDR3 rpm (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_cdr3), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA CDR3 reads (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_clonotypes), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA CDR3 clonotypes (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)
ggplot(rna_pos_results, aes(log10(tcrb_reads), log10(rna_clonotypes/rna_rpm), colour=cohort)) + geom_point() + labs(x="TCRB RNA reads (log10)", y="RNA clonal diversity (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

exome_rna_pos_results <- filter(results, exome_clonotypes>0, rna_clonotypes>0)
ggplot(exome_rna_pos_results, aes(log10(exome_clonotypes/exome_rpm), log10(rna_clonotypes/rna_rpm), colour=cohort)) + geom_point() + labs(x="Exome clonal diversity (log10)", y="RNA clonal diversity (log10)") + theme(text=element_text(size=16)) + scale_color_viridis(discrete=TRUE)

cor.test(exome_pos_results$exome_reads, exome_pos_results$exome_clonotypes)
cor.test(rna_pos_results$tcrb_rpm, rna_pos_results$rna_rpm)
cor.test(rna_pos_results$tcrb_reads, rna_pos_results$rna_clonotypes)