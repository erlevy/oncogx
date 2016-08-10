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

# contrast infiltration of cancer types
exome_split <- lapply(cancer_split, '[[', "exome_rpm")
rna_split <- lapply(cancer_split, '[[', "rna_rpm")
blood_split <- lapply(cancer_split, '[[', "blood_rpm")
idna_split <- lapply(cancer_split, '[[', "iDNA_score")

#ggplot(cancer_current, aes(x=exome_rpm, y=rna_rpm, colour=cohort)) + geom_point()

#exome_avg <- sapply(exome_split, mean, na.rm=TRUE)
#rna_avg <- sapply(rna_split, mean, na.rm=TRUE)
exome_avg <- log(sapply(exome_split, mean, na.rm=TRUE), 10)
rna_avg <- log(sapply(rna_split, mean, na.rm=TRUE), 10)
cohort_sizes <- table(cancer_current$cohort)
avgs_plot <- as.data.frame(cbind(exome_avg, rna_avg, cohort_sizes))

### FIGURE: average RNA vs. tumor exome RPM by cohort ###

#ggplot(avgs_plot, aes(exome_avg, rna_avg, size=cohort_sizes, colour=rownames(avgs_plot))) + geom_point() + labs(title="Cohort average exome and RNA RPM", x="log(exome CDR3 rpm)", y="log(RNA CDR3 rpm)")
ggplot(avgs_plot, aes(exome_avg, rna_avg)) + geom_point() + labs(x="exome CDR3 rpm (log10)", y="RNA CDR3 rpm (log10)") + geom_text_repel(data=avgs_plot, aes(label=rownames(avgs_plot)), size=6) + theme(text=element_text(size=16))

### FIGURE: Venn diagram/overall profiling of cohort sizes ###

# overall profiling
exome_n <- nrow(filter(cancer_current, !is.na(exome_rpm)))
rna_n <- nrow(filter(cancer_current, !is.na(rna_rpm)))
blood_n <- nrow(filter(cancer_current, !is.na(blood_rpm)))
exome_rna_n <- nrow(filter(cancer_current, !is.na(rna_rpm) & !is.na(exome_rpm)))
exome_blood_n <- nrow(filter(cancer_current, !is.na(blood_rpm) & !is.na(exome_rpm)))
blood_rna_n <- nrow(filter(cancer_current, !is.na(rna_rpm) & !is.na(blood_rpm)))
exome_rna_blood_n <- nrow(filter(cancer_current, !is.na(rna_rpm) & !is.na(exome_rpm) & !is.na(blood_rpm)))

grid.newpage()
draw.triple.venn(area1 = exome_n, area2 = rna_n, area3 = blood_n, 
                n12 = exome_rna_n, n23 = blood_rna_n, n13 = exome_blood_n,
                n123 = exome_rna_blood_n, fill = c("skyblue", "orange", "green"),
                category = c("Tumor exome", "Tumor RNA", "Blood exome"), cex=2, cat.cex=2, 
                cat.pos=c(-10,10,180), cat.dist=c(0.05, 0.05, 0.05))

exome_rna_blood <- filter(cancer_current, !is.na(rna_rpm) & !is.na(exome_rpm) & !is.na(blood_rpm))

# initial numbers
exome_pos <- filter(cancer_current, exome_rpm>0)
rna_pos <- filter(cancer_current, rna_rpm>0)
blood_pos <- filter(cancer_current, blood_rpm>0)

### FIGURE: cohort size barplot ###

cohort_sizes <- as.data.frame(table(exome_rna_blood$cohort))
ggplot(cohort_sizes, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity") + labs(x="Cohort", y="Number of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

### FIGURE: compare read counts of tumor exome and RNA by cohort ###

# contrast read counts of cancer types
exome_read_split <- lapply(cancer_split, '[[', "exome_reads")
rna_read_split <- lapply(cancer_split, '[[', "rna_reads")
percent_tils_split <- lapply(cancer_split, '[[', "percent_tils")
tcrb_rpm_split <- lapply(cancer_split, '[[', "tcrb_rpm")
tcrb_rpm_log_split <- lapply(tcrb_rpm_split, FUN=log)
gsva_cluster_split <- lapply(cancer_split, '[[', "gsva_cluster")

exome_read_melt <- melt(exome_read_split)
exome_read_melt$value <- log(exome_read_melt$value, 10)
ggplot(exome_read_melt, aes(x=L1, y=value)) + geom_boxplot() + labs(title="Exome read counts by cohort", x="Cohort", y="log(read count)")

rna_read_melt <- melt(rna_read_split)
rna_read_melt$value <- log(rna_read_melt$value, 10)
ggplot(rna_read_melt, aes(x=L1, y=value)) + geom_boxplot() + labs(title="RNA read counts by cohort", x="Cohort", y="log(read count)")

# iDNA vs purity for cancer types
frac_pure<-aggregate(patient_purity~iDNA_score,data=cancer_current,FUN=function(x) length(x[x>=0.80])/length(x))
ggplot(frac_pure,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),patient_purity))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥80% purity")+theme_classic(base_size = 24)

# exome vs. blood
# run clonotypes_analysis.R first
exome_only <- filter(clonotypes, exome > 0)
exome_pos_n <- length(unique(exome_pos$aa))
blood_pos_n <- length(unique(blood_pos$aa))
both_n <- length(unique(intersect(exome_pos$aa, blood_pos$aa)))

# distribution of cohorts with "sticky" clonotype
shared_patients <- exome_blood_shared_patients_only_c66
shared_patients_results <- filter(cancer_current, patient_uuid %in% shared_patients)
exome_blood_or_n_patients <- filter(cancer_current, !is.na(exome_rpm) | !is.na(blood_rpm))$patient_uuid
exome_blood_baseline <- filter(cancer_current, patient_uuid %in% exome_blood_or_n_patients)
shared_patients_df <- data.frame(table(shared_patients_results$cohort), table(exome_blood_baseline$cohort))
shared_patients_df <- filter(shared_patients_df, Freq>0)
colnames(shared_patients_df)[1] <- "Cohort"
shared_patients_df$Cohort <- factor(shared_patients_df$Cohort)
shared_patients_df$Freq.2 <- shared_patients_df$Freq/shared_patients_df$Freq.1

### FIGURE: c66 cohort fraction barplot ###
ggplot(shared_patients_df, aes(x=reorder(Cohort, -Freq.2), y=Freq.2)) + geom_bar(stat="identity") + labs(x="Cohort", y="Fraction of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

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

# profile gender with c66
shared_patients_gender_df <- data.frame(table(shared_patients_results$gender), table(exome_blood_baseline$gender))
shared_patients_gender_table <- matrix(c(303, 4262, 45, 4485), nrow=2, byrow=TRUE)
non_shared_patients_gender_df <- data.frame(table(non_shared_patients_results$gender), table(exome_blood_baseline$gender))
non_shared_patients_gender_table <- matrix(c(124, 4441, 95, 4435), nrow=2, byrow=TRUE)

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

# exome CDR3 high without RNA
exome_pos_rna_neg <- filter(cancer_current,D=="D+",R=="R-")
exome_pos_rna_pos <- filter(cancer_current,D=="D+",R=="R+")
exome_rpm_list <- list(rna_neg=-log(exome_pos_rna_neg$exome_rpm,10), rna_pos=-log(exome_pos_rna_pos$exome_rpm))
exome_rpm_melt <- melt(exome_rpm_list)
#ggplot(exome_rpm_melt, aes(x=L1, y=value)) + geom_boxplot() + labs(title="Exome RPM by RNA CDR3 nonzero", x="RNA CDR3", y="-log(exome rpm)")

exome_med <- median(filter(cancer_current, exome_rpm>0)$exome_rpm)
exome_pos_rna_neg_high <- filter(exome_pos_rna_neg, exome_rpm > exome_med)
exome_pos_rna_neg_df <- data.frame(table(exome_pos_rna_neg_high$cohort), table(exome_pos_rna_neg$cohort))
colnames(exome_pos_rna_neg_df)[1] <- "Cohort"
#ggplot(clonotypes_cohort_df, aes(x=Var1, y=Freq)) + geom_bar(stat="identity") + labs(title="Percent of cohort that has a public clonotype", x="Cohort", y="Percent")
ggplot(exome_pos_rna_neg_df, aes(Freq, Freq.1)) + geom_point() + labs(title="Number of D+R- patients with any exome RPM vs high by cohort", x="Number of patients with high exome CDR3", y="Number of D+R- patients") + geom_text_repel(data=exome_pos_rna_neg_df, aes(label=Cohort))

exome_pos_rna_neg_frac_df <- data.frame(table(exome_pos_rna_neg_high$cohort)/table(exome_pos_rna_neg$cohort))
colnames(exome_pos_rna_neg_df)[1] <- "Cohort"
ggplot(exome_pos_rna_neg_frac_df, aes(x=Var1, y=Freq)) + geom_bar(stat="identity") + labs(title="Fraction of D+R- patients with high exome RPM by cohort", x="Cohort", y="Percent")

exome_pos_vs_rna_df <- data.frame(table(exome_pos_rna_neg$cohort), table(exome_pos_rna_pos$cohort))
exome_pos_vs_rna_df <- filter(exome_pos_vs_rna_df, Freq>0)
exome_pos_vs_rna_df$Var1.1 <- NULL
colnames(exome_pos_vs_rna_df)[1] <- "Cohort"
colnames(exome_pos_vs_rna_df)[2] <- "D+R-"
colnames(exome_pos_vs_rna_df)[3] <- "D+R+"
exome_pos_vs_rna_df$Cohort <- factor(exome_pos_vs_rna_df$Cohort)
exome_pos_vs_rna_df_melt <- melt(exome_pos_vs_rna_df)

ggplot(exome_pos_vs_rna_df_melt, aes(x=Cohort, y=value, fill=variable)) + geom_bar(stat="identity", position="fill") + labs(x="Cohort", y="Fraction") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

### FIGURE: non c66 cohort fraction barplot ###
ggplot(non_shared_patients_df, aes(x=reorder(Cohort, -Freq.2), y=Freq.2)) + geom_bar(stat="identity") + labs(x="Cohort", y="Fraction of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

# rna vs. dna

# gsva table
gsva_table <- sapply(gsva_cluster_split, table)
gsva_percent <- prop.table(gsva_table, 2)
#barplot(gsva_percent) #, col=c("darkblue", "red"), legend=c("ifng low", "ifng high"), ylab="Percent high/low")
gsva_table_df <- melt(gsva_table)
ggplot(gsva_table_df, aes(x=Var2, y=value, fill=Var1)) + geom_bar(stat="identity", position="fill") + labs(title="Fraction of patients by GSVA immune signature groupings", x="Cohort", y="Fraction") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# gsva facet plots
cancer_heatmap <- cancer_current[!is.na(cancer_current$T.regs),]
cancer_heatmap$gsva_annotated <- as.character(cancer_heatmap$gsva_cluster)
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="1","High",as.character(gsva_annotated)))
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="2","Low",as.character(gsva_annotated)))
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="3","Mixed_Innate",as.character(gsva_annotated)))
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="4","Mixed_Adaptive",as.character(gsva_annotated)))

# to exclude: ACC, COAD, DLBC, KICH, READ, TGCT
# most interesting ones: BLCA, BRCA, CESC, HNSC, LUAD, LUSC, OV, PAAD, STAD
cancer_heatmap$gsva_annotated <- factor(cancer_heatmap$gsva_annotated, levels=c("Low", "Mixed_Innate", "Mixed_Adaptive", "High"))
cancer_heatmap <- filter(cancer_heatmap, cohort == "BLCA" | cohort == "BRCA" | cohort == "CESC" |
                         cohort == "HNSC" | cohort == "LUAD" | cohort == "LUSC" |
                         cohort == "OV" | cohort == "PAAD" | cohort == "STAD")

### FIGURE: gsva vs. cdr3 rpm facet plots ###
ggplot(cancer_heatmap, aes(x=gsva_annotated, y=log(exome_rpm,10))) + geom_boxplot() + labs(x="", y="Tumor exome CDR3 RPM (log10)") + facet_wrap(~cohort, ncol=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) # + theme(text=element_text(size=16))
ggplot(cancer_heatmap, aes(x=gsva_annotated, y=log(rna_rpm,10))) + geom_boxplot() + labs(x="", y="Tumor RNA CDR3 RPM (log10)") + facet_wrap(~cohort, ncol=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))# + theme(text=element_text(size=16))
#ggplot(cancer_heatmap, aes(x=gsva_annotated, y=log(blood_rpm,10))) + geom_boxplot() + labs(x="", y="Blood exome CDR3 RPM (log10)") + facet_wrap(~cohort, ncol=3) # + theme(text=element_text(size=16))


# gsva
cancer_heatmap <- cancer_current[!is.na(cancer_current$T.regs),]

signatures <- select(cancer_heatmap, cohort, T.regs:Neut)
cohorts <- signatures[,1]
signatures <- signatures[,-1]

signature_cols <- c("T Regs", "T CD8", "T CD4 Naive", "T CD4 Mem. Resting", "T CD4 Mem. Activated",
                    "T Folicular Helper", "T Gamma Delta", "B Naive", "B Memory", "Plasma", "NK Resting",
                    "NK Activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
                    "Dendritic Resting", "Dendritic Activated", "Mast Resting", "Mast Activated",
                    "Eosinophils", "Neutrophils")

# colors still not working, need to copy work form "brca_final_paper.R")
cohort_colors <- color.scale(as.numeric(cohorts), color.spec="bluerange")
gsva_colors <- color.scale(as.numeric(cancer_heatmap$gsva_cluster), color.spec="bluerange")
cohort_colors <- matrix(c(cohort_colors, cohort_colors), ncol=2, byrow=FALSE)

heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none", Colv=NA, labRow="",
             labCol=signature_cols, RowSideColors=cohort_colors)

#heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none", Colv=NA, labRow="",
#             labCol=signature_cols, cexCol=2, margins=c(22,6))

# clonal diversity


# survival
surv_results_exome <- function(surv_matrix)
{
  surv <- surv_matrix[,c("exome_level", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  #  plot(my.fit, col=c(1,2), main="Survival of HR+ by CDR3 groups", xlab="Survival Time (years)", ylab="Survival Probability")
  #  legend(0.5,0.4, c("nonzero", "zero"), fill=c(1,2), cex=0.75)
  coxph(my.surv ~ surv[,1], method="breslow")  
}

surv_results_RNA <- function(surv_matrix)
{
  surv <- surv_matrix[,c("RNA_level", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  #  plot(my.fit, col=c(1,2), main="Survival of HR+ by CDR3 groups", xlab="Survival Time (years)", ylab="Survival Probability")
  #  legend(0.5,0.4, c("nonzero", "zero"), fill=c(1,2), cex=0.75)
  coxph(my.surv ~ surv[,1], method="breslow")  
}

types <- split(cancer_current, cancer_current$cohort)
types_surv_exome <- lapply(types, surv_results_exome)
types_surv_RNA <- lapply(types, surv_results_RNA)

### FIGURE: COAD survival by exome CDR3 RPM ###
pdf("coad_exome_cdr3_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$COAD
surv <- surv_matrix[,c("exome_level", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="bottomright", c("Zero", "Nonzero"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.0006", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

### FIGURE: CESC survival by RNA CDR3 RPM ###
pdf("cesc_rna_cdr3_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$CESC
surv <- surv_matrix[,c("RNA_level", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="bottomright", c("Low", "High"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.004", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

### OLD 
# figure FracRNA pos vs DNA
frac_withRNA<-aggregate(rna_rpm~iDNA_score,data=filter(cdr3_results_current,!is.na(rna_rpm)),FUN=function(x) length(x[x>0.0255])/length(x))
ggplot(frac_withRNA,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors \n with high RNA CDR3 RPM")+theme_classic(base_size = 24)

#figure FracTIL vs DNA
frac_withTILs<-aggregate(percent_tils~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=5])/length(x))
ggplot(frac_withTILs,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),percent_tils))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥5% TILs")+theme_classic(base_size = 24)

# distribution by clinical groups
ggplot(cdr3_results_current,aes(as.factor(iDNA_score),fill=cohort))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())

# survival

