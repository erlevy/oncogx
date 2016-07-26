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

results <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-13-2016.txt", sep="\t")
clonotypes <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_clonotypes_all_7-14-2016.txt", sep="\t")

### recreate figures from BRCA paper

cdr3_results_current<-results

# calculating the iDNA score, rna_rpm, R and D groups
cdr3_results_current<-mutate(cdr3_results_current,exome_rpm=exome_cdr3*10^6/exome_reads)
tmp<-select(filter(cdr3_results_current,exome_rpm!=0),patient_uuid,exome_rpm)
tmp<-mutate(tmp,iDNA_score=ntile(exome_rpm,10))
cdr3_results_current<-left_join(cdr3_results_current,tmp)
cdr3_results_current<-mutate(cdr3_results_current,iDNA_score=ifelse(is.na(iDNA_score),0,iDNA_score))

cdr3_results_current<-mutate(cdr3_results_current,rna_rpm=rna_cdr3*10^6/rna_reads)
cdr3_results_current<-mutate(cdr3_results_current,blood_rpm=blood_cdr3*10^6/blood_reads)
cdr3_results_current<-mutate(cdr3_results_current,tcrb_rpm=tcrb_reads*10^6/rna_reads)
cdr3_results_current<-mutate(cdr3_results_current,D=ifelse(exome_rpm==0,"D-","D+"),R=ifelse(rna_rpm==0,"R-","R+"))

### recreate relevant results report
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

ggplot(avgs_plot, aes(exome_avg, rna_avg, size=cohort_sizes, colour=rownames(avgs_plot))) + geom_point() + labs(title="Cohort average exome and RNA RPM", x="log(exome CDR3 rpm)", y="log(RNA CDR3 rpm)")
ggplot(avgs_plot, aes(exome_avg, rna_avg)) + geom_point() + labs(title="Cohort average exome and RNA RPM", x="log(exome CDR3 rpm)", y="log(RNA CDR3 rpm)") + geom_text_repel(data=avgs_plot, aes(label=rownames(avgs_plot)))

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

# exome vs. blood
# run clonotypes_analysis.R first
exome_only <- filter(clonotypes, exome > 0)

exome_pos_n <- length(unique(exome_pos$aa))
blood_pos_n <- length(unique(blood_pos$aa))
both_n <- length(unique(intersect(exome_pos$aa, blood_pos$aa)))

public_clonotypes_results <- filter(cancer_current, patient_uuid %in% exome_blood_public_patients)
public_clonotypes_baseline <- filter(cancer_current, patient_uuid %in% exome_or_blood_patients)
table(public_clonotypes_results$cohort)

clonotypes_cohort_df <- data.frame(table(public_clonotypes_results$cohort), table(public_clonotypes_baseline$cohort))
colnames(clonotypes_cohort_df)[1] <- "Cohort"
#ggplot(clonotypes_cohort_df, aes(x=Var1, y=Freq)) + geom_bar(stat="identity") + labs(title="Percent of cohort that has a public clonotype", x="Cohort", y="Percent")
ggplot(clonotypes_cohort_df, aes(Freq, Freq.1)) + geom_point() + labs(title="Number of patients with any CDR3 reads vs. public clonotypes", x="Number of patients with public clonotypes", y="Number of patients with CDR3 reads") + geom_text_repel(data=clonotypes_cohort_df, aes(label=Cohort))

grid.newpage()
draw.pairwise.venn(area1 = exome_pos_n, area2 = blood_pos_n, cross.area = both_n, fill = c("skyblue", "orange"))

# exome CDR3 high without RNA
tmp <- filter(cancer_current,D=="D+",R=="R-")

# rna vs. dna

# gsva table
gsva_table <- sapply(gsva_cluster_split, table)
gsva_percent <- prop.table(gsva_table, 2)
barplot(gsva_percent) #, col=c("darkblue", "red"), legend=c("ifng low", "ifng high"), ylab="Percent high/low")

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

