library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)
library(car)
library(gvlma)
library(MASS)
library(leaps)
library(relaimpo)
library(dplyr)

### January 29, 2016 ###
# Purpose is to generate final tables and figures for the BRCA TCRSEQ paper
# Using finalized dataset that has 1078 samples and NAs for any missing fields

# Input: final results table to be used in paper
#results <- read.csv("cdr3_results_with_expression_blood_ptprc_1078.txt", sep="\t")
results <- read.csv("cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_clustering_age_dedup_1078.txt", sep="\t")
# Contains: patient and sample information, exome and rna imseq information, lymphocyte, clinical, groupings,
# and the GSVA scores for the LM22 signature

exome_freq <- as.numeric(results$exome_imseq)/as.numeric(results$exome_reads)
rna_freq <- as.numeric(results$rna_imseq)/as.numeric(results$rna_reads)

lymph_per <- as.numeric(results$lymphocyte_presence)
exome_with_lymph <- exome_freq[which(lymph_per==1)]
exome_without_lymph <- exome_freq[which(lymph_per==0)]
rna_with_lymph <- rna_freq[which(lymph_per==1)]
rna_without_lymph <- rna_freq[which(lymph_per==0)]

# p-value of assocation
exome_wilcox <- wilcox.test(exome_freq ~ lymph_per)
exome_pval <- exome_wilcox$p.value
rna_wilcox <- wilcox.test(rna_freq ~ lymph_per)
rna_pval <- rna_wilcox$p.value

# convert to RPM
exome_rpm <- as.numeric(results$exome_imseq)/as.numeric(results$exome_reads)*1000000
rna_rpm <- as.numeric(results$rna_imseq)/as.numeric(results$rna_reads)*1000000

# Figure 2a: overall distribution of cdr3 reads across samples
plot(sort(exome_rpm), xlab="Patients", ylab="Exome CDR3 RPM")

### Getting numbers for the paper ###
# number of patients with CDR3 reads
cdr3_positive <- which(results$exome_imseq>0) # 473
# number of patients with CDR3 read with 0 TILs
cdr3_no_til <- intersect(which(results$exome_imseq>0), which(results$lymphocyte_percent==0))  # 144
# number of patients with >10% TILs and 0 CDR3 reads
til_10_no_cdr3 <- intersect(which(results$exome_imseq==0), which(results$lymphocyte_percent>10))  #65
# number of patients with RNA-seq
rna_patients <- which(!is.na(results$rna_imseq))  # 1074
# number of patients with RNA-seq clonotypes
rna_cdr3 <- which(results$rna_imseq>0) # 906
# number of patients with cdr3 in exome and RNA
cdr3_exome_rna <- intersect(which(results$exome_imseq>0), which(results$rna_imseq>0)) # 435

### clonotype analysis ###
# get shared clonotypes in the blood samples
shared <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
imseq <- matrix(nrow=0, ncol=3)
colnames(imseq) <- c("nucleotide_sequence", "aa_sequence", "number_of_patients")
blood_shared <- rep(0, nrow(results))
imseq_dir <- '/Users/Eric/tcga/clonotypes/final/blood/'
has_shared_blood <- c()
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  imseq_file <- paste(results$blood_id[i], ".tsv", sep="")
  imseq_path <- paste(imseq_dir, imseq_file, sep="")
  if (!file.exists(imseq_path))
  {
    has_shared_blood <- c(has_shared_blood, 0)
    next
  }
  if(file.info(imseq_path)$size!=0)
  {
    imseq_reads <- as.matrix(read.csv(imseq_path, sep="\t", header=T))
    imseq_reads <- imseq_reads[!duplicated(imseq_reads[,10]),, drop=FALSE] 
    has_shared_blood <- c(has_shared_blood, length(intersect(shared, imseq_reads[,10]))) 
  }
  if (nrow(imseq_reads)>0)
  {
    for (j in 1:nrow(imseq_reads))
    {
      nucleotide <- imseq_reads[j,10]
      if (nucleotide==shared) {blood_shared[i] <- 1}
      aa <- imseq_reads[j,11]
      imseq_match <- which(imseq[,1]==nucleotide)
      if (length(imseq_match)==0)
      {
        imseq <- rbind(imseq, c(nucleotide, aa, "1"))
      }
      else
      {
        imseq[imseq_match,3] <- as.character(as.numeric(imseq[imseq_match,3])+1)
      }
    }
  }
}
imseq_blood <- imseq
#write.table(imseq, "repeated_clonotypes_normal_1078.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# get shared clonotype in the tumor samples
shared <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
imseq <- matrix(nrow=0, ncol=3)
colnames(imseq) <- c("nucleotide_sequence", "aa_sequence", "number_of_patients")
tumor_shared <- rep(0, nrow(results))
imseq_dir <- '/Users/Eric/tcga/clonotypes/final/exome/'
has_shared_tumor <- c()
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  imseq_file <- paste(results$exome_id[i], ".tsv", sep="")
  imseq_path <- paste(imseq_dir, imseq_file, sep="")
  if (!file.exists(imseq_path))
  {
    next
  }
  if(file.info(imseq_path)$size!=0)
  {
    imseq_reads <- as.matrix(read.csv(imseq_path, sep="\t", header=T))
    imseq_reads <- imseq_reads[!duplicated(imseq_reads[,10]),, drop=FALSE] 
    has_shared_tumor <- c(has_shared_tumor, length(intersect(shared, imseq_reads[,10]))) 
  }
  if (nrow(imseq_reads)>0)
  {
    for (j in 1:nrow(imseq_reads))
    {
      nucleotide <- imseq_reads[j,10]
      if (nucleotide==shared) {tumor_shared[i] <- 1}
      aa <- imseq_reads[j,11]
      imseq_match <- which(imseq[,1]==nucleotide)
      if (length(imseq_match)==0)
      {
        imseq <- rbind(imseq, c(nucleotide, aa, "1"))
      }
      else
      {
        imseq[imseq_match,3] <- as.character(as.numeric(imseq[imseq_match,3])+1)
      }
    }
  }
}
imseq_tumor <- imseq
#write.table(imseq, "repeated_clonotypes_tumor_1078.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# get shared clonotypes in the rna samples
shared <- "TGTGCCACCAGCAGAGACACAGAGCTGCAGTGCTTCCTGCTCTCTGTTCATAAACCTCATTGTTTCCCAGATCCAGGTGCTTTCTCT"
imseq <- matrix(nrow=0, ncol=3)
colnames(imseq) <- c("nucleotide_sequence", "aa_sequence", "number_of_patients")
rna_shared <- rep(0, nrow(results))
imseq_dir <- '/Users/Eric/tcga/clonotypes/final/rna/'
tumor_dir <- '/Users/Eric/tcga/clonotypes/final/exome/'
rna_exome_shared <- c()
for (i in 1:nrow(results))
{
  imseq_reads <- 0
  imseq_file <- paste(results$rna_id[i], ".tsv", sep="")
  imseq_path <- paste(imseq_dir, imseq_file, sep="")
  tumor_file <- paste(results$exome_id[i], ".tsv", sep="")
  tumor_path <- paste(tumor_dir, tumor_file, sep="")
  if (!file.exists(imseq_path))
  {
    next
  }
  if(file.info(imseq_path)$size!=0)
  {
    imseq_reads <- as.matrix(read.csv(imseq_path, sep="\t", header=T))
    imseq_reads <- imseq_reads[!duplicated(imseq_reads[,10]),, drop=FALSE] 
    tumor_reads <- as.matrix(read.csv(tumor_path, sep="\t", header=T))
    tumor_reads <- tumor_reads[!duplicated(tumor_reads[,10]),, drop=FALSE] 
    rna_exome_shared <- c(rna_exome_shared, length(intersect(imseq_reads[,10], tumor_reads[,10])))
  }
  if (nrow(imseq_reads)>0)
  {
    for (j in 1:nrow(imseq_reads))
    {
      nucleotide <- imseq_reads[j,10]
      if (nucleotide==shared) {rna_shared[i] <- 1}
      aa <- imseq_reads[j,11]
      imseq_match <- which(imseq[,1]==nucleotide)
      if (length(imseq_match)==0)
      {
        imseq <- rbind(imseq, c(nucleotide, aa, "1"))
      }
      else
      {
        imseq[imseq_match,3] <- as.character(as.numeric(imseq[imseq_match,3])+1)
      }
    }
  }
}
imseq_rna <- imseq
#write.table(imseq, "repeated_clonotypes_rna_1078.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# clonotypes per read
patients_cdr3_positive <- results[which(results$exome_imseq>0),]
summary(patients_cdr3_positive$exome_clonotypes) # mean: 1.913, max: 13
# clonotypes in common between all of rna and exome
clonotypes_shared_all <- intersect(imseq_rna[,1], imseq_tumor[,1])  # 15
clonotypes_all <- union(imseq_rna[,1], imseq_tumor[,1])  # 7954

# patients with rna and exome shared reads
sum(rna_exome_shared) # 13
# patients with the shared clonotype in both blood and tumor
length(intersect(which(has_shared_blood==1), which(has_shared_tumor==1)))

# heatmap of gsva clusters
results_clustering <- results[!is.na(results$gsva_cluster),]
#results_clustering <- results_clustering[!is.na(results_clustering$rna_reads),]
#results_clustering <- results_clustering[which(results_clustering$clinical_group!="other"),]
#results_clustering <- results_clustering[!is.na(results_clustering$clinical_group),]
col_mat <- matrix(nrow=nrow(results_clustering), ncol=4)
rna_freq <- results_clustering$rna_imseq/results_clustering$rna_reads
rna_med <- median(rna_freq)
#for (i in 1:nrow(col_mat))
#{
  group <- as.character(results_clustering$clinical_group[i])
  exome <- results_clustering$exome_imseq[i]
  rna <- rna_freq[i]
  cluster_patient <- results_clustering$gsva_cluster[i]
  if (exome==0) {col_mat[i,2]="white"}
  else {col_mat[i,2]="black"}
  if (rna<=rna_med) {col_mat[i,3]="white"}
  else {col_mat[i,3]="purple"}
  if (group=="her2+") {col_mat[i,1]="red"}
  else if (group=="her2-") {col_mat[i,1]="green"}
  else if (group=="triple-") {col_mat[i,1]="blue"}
  if (cluster_patient==1) {col_mat[i,4]="red"}
  else if (cluster_patient==2) {col_mat[i,4]="green"}
  else if (cluster_patient==3) {col_mat[i,4]="blue"}
}

ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

for (i in 1:nrow(col_mat))
{
  exome <- results_clustering$exome_imseq[i]
  cluster_patient <- results_clustering$gsva_cluster[i]
#  cluster_patient <- results_clustering$mycl[i]
  if (exome==0) {col_mat[i,2]="white"}
  else {col_mat[i,2]="black"}
  if (cluster_patient==1) {col_mat[i,4]="#619CFF"}
  else if (cluster_patient==2) {col_mat[i,4]="#F8766D"}
  else if (cluster_patient==3) {col_mat[i,4]="#00BA38"}
  else if (cluster_patient==4) {col_mat[i,4]="grey"}
}

col_mat[,2] <- col_mat[,4]
col_mat[,3] <- col_mat[,4]
col_mat[1072,1] <- NA

signatures <- results_clustering[,37:58] # new cdr3_results with dedup
rownames(signatures) <- results_clustering[,1]
signatures <- scale(signatures)

signature_cols <- c("T Regs", "T CD8", "T CD4 Naive", "T CD4 Mem. Resting", "T CD4 Mem. Activated",
                    "T Folicular Helper", "T Gamma Delta", "B Naive", "B Memory", "Plasma", "NK Resting",
                    "NK Activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
                    "Dendritic Resting", "Dendritic Activated", "Mast Resting", "Mast Activated",
                    "Eosinophils", "Neutrophils")

heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none",
             Colv=NA, labRow="", cexCol=2.25, margins=c(30,8), labCol=signature_cols)

exome_rpm <- (results$exome_imseq/results$exome_reads)*1000000
rna_rpm <- (results$rna_imseq/results$rna_reads)*1000000
blood_rpm <- (results$blood_cdr3/results$blood_reads)*1000000
tcrb_rpm <- (results$tcrb_reads/results$rna_reads)*1000000 
results <- cbind(results, exome_rpm, rna_rpm, blood_rpm, tcrb_rpm)

# chi-squared test
cdr3_nonzero <- results$exome_imseq
cdr3_nonzero[which(cdr3_nonzero>0)] <- 1
cdr3_nonzero[which(cdr3_nonzero==0)] <- 0
cdr3_nonzero <- as.factor(cdr3_nonzero)
results <- cbind(results, cdr3_nonzero)
clusters_methy <- as.factor(results$CLUS_Methlyation_cNMF)
clusters_mrna <- as.factor(results$CLUS_mRNAseq_cNMF)
clusters_gsva <- as.factor(results$gsva_cluster)

chi_methy <- chisq.test(cdr3_nonzero, clusters_methy)
chi_mrna <- chisq.test(cdr3_nonzero, clusters_mrna)
chi_gsva <- chisq.test(cdr3_nonzero, clusters_gsva)

# fishers with clinical groups and cdr3/cytolytic
restriction_cdr3_columns <- c("cdr3_nonzero", "foxp3", "rate_non", "rna_rpm", "blood_rpm", "clinical_group")
cdr3_results_cdr3_restricted <- results[,restriction_cdr3_columns]
cdr3_results_cdr3_restricted <- cdr3_results_cdr3_restricted[complete.cases(cdr3_results_cdr3_restricted),]
cdr3_results_cdr3_restricted <- filter(cdr3_results_cdr3_restricted,clinical_group!="other")

cytolytic_group <- results$cytolytic
cytolytic_med <- median(cytolytic_group, na.rm=TRUE)
cytolytic_group[which(cytolytic_group<cytolytic_med)] <- 0
cytolytic_group[which(cytolytic_group>=cytolytic_med)] <- 1
cytolytic_group <- as.factor(cytolytic_group)
results <- cbind(results, cytolytic_group)
restriction_cytolytic_columns <- c("cytolytic_group", "tcrb_rpm", "cd4", "cd8a", "foxp3", "rna_rpm", "clinical_group")
cdr3_results_cytolytic_restricted <- results[,restriction_cytolytic_columns]
cdr3_results_cytolytic_restricted <- cdr3_results_cytolytic_restricted[complete.cases(cdr3_results_cytolytic_restricted),]
cdr3_results_cytolytic_restricted <- filter(cdr3_results_cytolytic_restricted,clinical_group!="other")

chi_clinical_cdr3 <- chisq.test(cdr3_results_cdr3_restricted$cdr3_nonzero, cdr3_results_cdr3_restricted$clinical_group)
chi_clinical_cytolytic <- chisq.test(cdr3_results_cytolytic_restricted$cytolytic_group, cdr3_results_cytolytic_restricted$clinical_group)

restriction_cytolytic_nonclinical_columns <- c("cytolytic", "tcrb_rpm", "cd4", "cd8a", "foxp3", "rna_rpm")
cdr3_results_cytolytic_nonclinical_restricted <- results[,restriction_cytolytic_nonclinical_columns]
cdr3_results_cytolytic_nonclinical_restricted <- cdr3_results_cytolytic_nonclinical_restricted[complete.cases(cdr3_results_cytolytic_nonclinical_restricted),]

# lm
lm_vars <- c("exome_rpm", "tcrb_rpm", "cd4", "cd8a", "foxp3", "cytolytic", "CLUS_mRNAseq_cNMF",
             "CLUS_Methlyation_cNMF", "rate_non", "patient_purity", "gsva_cluster", "rna_rpm",
             "blood_rpm", "clinical_group", "CLI_years_to_birth")
lm_nonfactor <- c("exome_rpm", "tcrb_rpm", "cd4", "cd8a", "foxp3", "cytolytic", "rate_non", 
                  "patient_purity", "rna_rpm", "blood_rpm","CLI_years_to_birth")
results_lm <- results[,lm_vars]
results_lm_nonfactor <- results[,lm_nonfactor]
#results_lm <- scale(results_lm)
results_lm <- results_lm[complete.cases(results_lm),]
results_lm_nonfactor <- results_lm_nonfactor[complete.cases(results_lm_nonfactor),]

rpm_lm_all <- lm(exome_rpm~tcrb_rpm+cd4+cd8a+foxp3+cytolytic+CLUS_mRNAseq_cNMF+CLUS_Methlyation_cNMF+
               rate_non+patient_purity+gsva_cluster+rna_rpm+blood_rpm+clinical_group+CLI_years_to_birth,
               data=results_lm)

rpm_lm <- lm(exome_rpm~tcrb_rpm+cd4+cd8a+foxp3+cytolytic+
               rate_non+patient_purity+gsva_cluster+rna_rpm+blood_rpm+clinical_group, data=results_lm)

rpm_lm <- lm(exome_rpm~foxp3+rate_non+rna_rpm+blood_rpm+clinical_group, data=results)

rpm_lm_nonfactor <- lm(exome_rpm~tcrb_rpm+cd4+cd8a+foxp3+cytolytic+rate_non+patient_purity+
                       rna_rpm+blood_rpm+CLI_years_to_birth,
                       data=results)

rpm_lm_nonfactor <- lm(exome_rpm~foxp3+rate_non+rna_rpm+blood_rpm, data=results)

# logit
nonzero_logit <- glm(cdr3_nonzero~tcrb_rpm+cd4+cd8a+foxp3+cytolytic+CLUS_mRNAseq_cNMF
                     +CLUS_Methlyation_cNMF+rate_non+patient_purity+gsva_cluster+rna_rpm+blood_rpm+
                       clinical_group+CLI_years_to_birth,
                     data=results, family='binomial')

nonzero_logit <- glm(cdr3_nonzero~tcrb_rpm+cd4+cd8a+foxp3+cytolytic+CLUS_mRNAseq_cNMF
                     +rate_non+patient_purity+gsva_cluster+rna_rpm+blood_rpm,
                     data=results, family='binomial')

nonzero_logit <- glm(cdr3_nonzero~foxp3+gsva_cluster+blood_rpm,
                     data=results, family='binomial')

# rna rpm as function of exome and tcrb
exome_rpm_c <- scale(exome_rpm, scale=FALSE)
tcrb_rpm_c <- scale(tcrb_rpm, scale=FALSE)
rna_lm <- lm(rna_rpm~exome_rpm_c*tcrb_rpm_c)

# variable selection
step_fit <- stepAIC(rpm_lm_all, direction="both")
step_fit$anova

# exome_rpm ~ foxp3 + rate_non + rna_rpm + blood_rpm + clinical_group

#leaps <- regsubsets(exome_rpm~tcrb_rpm+cd4+cd8a+foxp3+cytolytic+CLUS_mRNAseq_cNMF+CLUS_Methlyation_cNMF+
#                    rate_non+patient_purity+gsva_cluster+rna_rpm+blood_rpm+clinical_group+CLI_years_to_birth,
#                    data=results_lm, nbest=10)
#plot(leaps, scale="r2")
#subsets(leaps, statistic="rsq")

# relative importance
calc <- calc.relimp(rpm_lm, type = c("lmg", "last", "first"), rela = FALSE)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(rpm_lm_nonfactor, b = 1000, type = c("lmg", "last", "first", "pratt"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)

booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result

# testing cytolytic activity
cyto_lm_all <- lm(cytolytic~exome_rpm+tcrb_rpm+cd4+cd8a+foxp3+CLUS_mRNAseq_cNMF+CLUS_Methlyation_cNMF+
                 rate_non+patient_purity+gsva_cluster+rna_rpm+blood_rpm+clinical_group+CLI_years_to_birth,
                 data=results_lm)

# variable selection
step_fit <- stepAIC(cyto_lm_all, direction="both")
step_fit$anova

# final models
cyto_lm_oh <- lm(cytolytic~exome_rpm+tcrb_rpm+cd8a+clinical_group+CLI_years_to_birth,
                 data=results)

cyto_lm_aic <- lm(cytolytic~tcrb_rpm+cd4+cd8a+foxp3+rna_rpm+clinical_group,
                 data=results)
