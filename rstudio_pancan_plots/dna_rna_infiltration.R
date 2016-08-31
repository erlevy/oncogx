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
library(FSA)

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
#cancer_current <- filter(cdr3_results_current, cohort != "THYM", cohort != "UVM", cohort != "CHOL")
cancer_current <- filter(cdr3_results_current, cohort != "THYM", cohort != "UVM", cohort != "CHOL",
                         cohort != "ACC", cohort != "CESC", cohort != "DLBC", cohort != "ESCA", 
                         cohort != "KICH", cohort != "LGG", cohort != "LIHC", cohort != "MESO",
                         cohort != "PCPG", cohort != "TGCT", cohort != "UCS")
cancer_current$cohort <- factor(cancer_current$cohort)
cancer_split <- split(cancer_current, cancer_current$cohort)

# medians of nonzero pancan
nonzero_pancan_exome <- filter(cancer_current, exome_cdr3 > 0)
nonzero_pancan_rna <- filter(cancer_current, rna_cdr3 > 0)
nonzero_pancan_blood <- filter(cancer_current, blood_cdr3 > 0)
nonzero_pancan_tils <- filter(cancer_current, percent_tils > 0)

median(nonzero_pancan_exome$exome_cdr3)
quantile(nonzero_pancan_exome$exome_cdr3, 0.9)
max(nonzero_pancan_exome$exome_cdr3)

median(nonzero_pancan_rna$rna_cdr3)
quantile(nonzero_pancan_rna$rna_cdr3, 0.9)
max(nonzero_pancan_rna$rna_cdr3)

median(nonzero_pancan_blood$blood_cdr3)
quantile(nonzero_pancan_blood$blood_cdr3, 0.9)
max(nonzero_pancan_blood$blood_cdr3)

median(nonzero_pancan_tils$percent_tils)
quantile(nonzero_pancan_tils$percent_tils, 0.9)
max(nonzero_pancan_tils$percent_tils)

nrow(nonzero_pancan_exome)
nrow(filter(cancer_current, exome_cdr3==0))

nrow(nonzero_pancan_rna)
nrow(filter(cancer_current, rna_cdr3==0))

nrow(nonzero_pancan_blood)
nrow(filter(cancer_current, blood_cdr3==0))

nrow(nonzero_pancan_tils)
nrow(filter(cancer_current, percent_tils==0))

# correlations with TILs and purity
cor.test(cancer_current$exome_rpm, cancer_current$percent_tils, method="spearman")
cor.test(cancer_current$rna_rpm, cancer_current$percent_tils, method="spearman")
cor.test(cancer_current$blood_rpm, cancer_current$percent_tils, method="spearman")

cor.test(cancer_current$exome_rpm, cancer_current$patient_purity, method="spearman")
cor.test(cancer_current$rna_rpm, cancer_current$patient_purity, method="spearman")
cor.test(cancer_current$blood_rpm, cancer_current$patient_purity, method="spearman")

### non-parametric exome testing
# kruskal-wallis by cohort
krusk <- kruskal.test(exome_rpm ~ cohort, data=cancer_current)
dunn <- dunnTest(exome_rpm ~ cohort, data=cancer_current, method="bh")

# create dunn matrix
dunn_mat <- matrix(NA, nrow=length(unique(cancer_current$cohort)), ncol=length(unique(cancer_current$cohort)))
rownames(dunn_mat) <- sort(unique(cancer_current$cohort))
colnames(dunn_mat) <- sort(unique(cancer_current$cohort))
for (i in 1:nrow(dunn_mat)) {
  row_name <- rownames(dunn_mat)[i]
  for (j in 1:i) {
    col_name <- colnames(dunn_mat)[j]
    if (row_name == col_name) {dunn_mat[i,j] <- 1}
    else {
      dun_access <- paste(col_name, "-")
      dun_access <- paste(dun_access, row_name)
      dun_i <- which(as.character(dunn$res$Comparison)==dun_access) 
      dun_p <- dunn$res$P.adj[dun_i]
      dunn_mat[i,j] <- dun_p
      dunn_mat[j,i] <- dun_p
    }
  }
}

heatmap.plus(-log10(dunn_mat), col=bluered(51), symm=TRUE, keep.dendro=FALSE, Rowv=NA, Colv=NA)

arrange(dunn$res, P.adj)


### non-parametric RNA testing
# kruskal-wallis by cohort
krusk <- kruskal.test(rna_rpm ~ cohort, data=cancer_current)
dunn <- dunnTest(rna_rpm ~ cohort, data=cancer_current, method="bh")

# create dunn matrix
dunn_mat <- matrix(NA, nrow=length(unique(cancer_current$cohort)), ncol=length(unique(cancer_current$cohort)))
rownames(dunn_mat) <- sort(unique(cancer_current$cohort))
colnames(dunn_mat) <- sort(unique(cancer_current$cohort))
for (i in 1:nrow(dunn_mat)) {
  row_name <- rownames(dunn_mat)[i]
  for (j in 1:i) {
    col_name <- colnames(dunn_mat)[j]
    if (row_name == col_name) {dunn_mat[i,j] <- 1}
    else {
      dun_access <- paste(col_name, "-")
      dun_access <- paste(dun_access, row_name)
      dun_i <- which(as.character(dunn$res$Comparison)==dun_access) 
      dun_p <- dunn$res$P.adj[dun_i]
      dunn_mat[i,j] <- dun_p
      dunn_mat[j,i] <- dun_p
    }
  }
}

heatmap.plus(-log10(dunn_mat), col=bluered(51), symm=TRUE, keep.dendro=FALSE, Rowv=NA, Colv=NA)

arrange(dunn$res, P.adj)

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

exome_not_na <- filter(cancer_current, !is.na(exome_rpm))
rna_not_na <- filter(cancer_current, !is.na(rna_rpm))
blood_not_na <- filter(cancer_current, !is.na(blood_rpm))

### FIGURE: cohort size barplot ###

cohort_sizes <- as.data.frame(table(exome_rna_blood$cohort))
ggplot(cohort_sizes, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity") + labs(x="Cohort", y="Number of patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=16))

### FIGURE: compare read counts of tumor exome and RNA by cohort ###

# contrast read counts of cancer types
exome_read_split <- lapply(cancer_split, '[[', "exome_reads")
exome_rpm_split <- lapply(cancer_split, '[[', "exome_rpm")
rna_read_split <- lapply(cancer_split, '[[', "rna_reads")
rna_rpm_split <- lapply(cancer_split, '[[', "rna_rpm")
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

exome_rpm_melt <- melt(exome_rpm_split)
exome_rpm_melt$value <- log(exome_rpm_melt$value, 10)
ggplot(exome_rpm_melt, aes(x=L1, y=value)) + geom_boxplot() + labs(title="Tumor exome RPM counts by cohort", x="Cohort", y="log(read count)")

rna_rpm_melt <- melt(rna_rpm_split)
rna_rpm_melt$value <- log(rna_rpm_melt$value, 10)
ggplot(rna_rpm_melt, aes(x=L1, y=value)) + geom_boxplot() + labs(title="Tumor rna RPM counts by cohort", x="Cohort", y="log(read count)")

# iDNA vs purity for cancer types
frac_pure<-aggregate(patient_purity~iDNA_score,data=cancer_current,FUN=function(x) length(x[x>=0.80])/length(x))
ggplot(frac_pure,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),patient_purity))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥80% purity")+theme_classic(base_size = 24)


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

# check D+ for tcrb vs. R+/R-
exome_pos_new <- mutate(exome_pos,R2=ifelse(rna_rpm==0,0,1))

wilcox.test(exome_pos_new$tcrb_rpm, exome_pos_new$R2)
boxplot(lapply(split(exome_pos$tcrb_rpm, exome_pos$R), log10))
rna_tcrb_box <- melt(lapply(split(exome_pos$tcrb_rpm, exome_pos$R), log10))
ggplot(rna_tcrb_box, aes(x=L1, y=value)) + geom_boxplot() + labs(x="RNA CDR3", y="TCRB RPM (log10)") + theme(text=element_text(size=16))

rna_exome_box <- melt(lapply(split(exome_pos$exome_rpm, exome_pos$R), log10))
ggplot(rna_exome_box, aes(x=L1, y=value)) + geom_boxplot() + labs(x="RNA CDR3", y="Tumor Exome RPM (log10)") + theme(text=element_text(size=16))

# idna by cohort
ggplot(cancer_current,aes(as.factor(iDNA_score)))+geom_bar()+facet_wrap(~cohort, ncol=4, scales="free_y")

# exome high by cohort
exome_rpm_med <- median(filter(cancer_current, exome_level>0)$exome_rpm)
cancer_current<-mutate(cancer_current,exome_lowhigh=ifelse(exome_rpm>exome_rpm_med,"1","0"))
exome_nonzero <- filter(cancer_current, exome_level==1)

ggplot(cancer_current,aes(as.factor(iDNA_score)))+geom_bar()+facet_wrap(~cohort, ncol=4, scales="free_y")
ggplot(exome_nonzero,aes(as.factor(exome_lowhigh)))+geom_bar()+facet_wrap(~cohort, ncol=4, scales="free_y")


# gender comparison by cohort
compare_gender <- function(df_in)
{
  gender_split <- split(df_in, df_in$gender)  
  gender_fc <- mean(gender_split$male$rna_rpm, na.rm=TRUE)/mean(gender_split$female$rna_rpm, na.rm=TRUE)
}

gender_comparisons <- lapply(cancer_split, compare_gender)

