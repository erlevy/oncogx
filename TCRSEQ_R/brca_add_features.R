library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)
library(ggplot2)
library(dplyr)

results <- read.csv("cdr3_results_with_expression_blood_ptprc_1078.txt", sep="\t")

### expression
expression <- read.csv("BRCA_expression_firehose_patient_ids.csv")

# process expression matrix
expression <- as.matrix(expression)
rownames(expression) <- expression[,1]
expression <- expression[,-1]
mode(expression) <- 'numeric'
expression <- t(expression)
rownames(expression) <- gsub("X", "", fixed=TRUE, rownames(expression))
rownames(expression) <- gsub(".", "-", fixed=TRUE, rownames(expression))

# CD4, CD8A, CD20, FOXP3
# CD20 not there
patient_expression <- c()
gene_list <- c("CD4", "CD8A", "FOXP3", "GZMA", "PRF1")
for (i in 1:nrow(results))
{
  patient <- as.character(results$patient_id[i])
  matching <- which(rownames(expression)==patient)
  if (length(matching)>0)
  {
    matching_expression <- expression[matching, gene_list]
  }
  else
  {
    matching_expression <- rep(NA, length(gene_list))
  }
  patient_expression <- rbind(patient_expression, matching_expression)
}
rownames(patient_expression) <- seq(1,nrow(patient_expression))
results_expression <- cbind(results, patient_expression)
cytolytic <- sqrt(results_expression$GZMA*results_expression$PRF1)
results_expression <- cbind(results_expression, cytolytic)
results_all <- results_expression

### mutation
mutation <- read.csv("BRCA-TP.samplefeatures.txt", sep="\t")

# process mutation matrix
mutation <- as.matrix(mutation)
rownames(mutation) <- mutation[,1]
mutation <- mutation[,-1]

# match samples
patient_mutation <- c()
data_list <- c("CLUS_mRNA_cNMF", "CLUS_mRNAseq_cNMF", "CLUS_Methlyation_cNMF", 
               "rate_non", "rate_sil", "CLI_years_to_birth")
for (i in 1:nrow(results))
{
  sample_id <- substr(as.character(results$sample_id[i]),1,12)
  matching <- which(rownames(mutation)==sample_id)
  if (length(matching)>0)
  {
    matching_mutation <- mutation[matching, data_list]
  }
  else
  {
    matching_mutation <- rep(NA, length(data_list))
  }
  patient_mutation <- rbind(patient_mutation, matching_mutation)
}
rownames(patient_mutation) <- seq(1,nrow(patient_mutation))
mode(patient_mutation) <- 'numeric'
results_mutation <- cbind(results, patient_mutation)
results_clustering <- results_mutation
results_all <- cbind(results_all, patient_mutation)

### purity
purity <- read.csv("/Users/Eric/BRCA/pancancer_purity_butte.txt", sep="\t")

# process purity matrix
purity <- as.matrix(purity)
rownames(purity) <- purity[,1]
purity <- purity[,-1]

# match samples
patient_purity <- c()
data_list <- c("CPE")
for (i in 1:nrow(results))
{
  sample_id <- substr(as.character(results$sample_id[i]),1,16)
  matching <- which(rownames(purity)==sample_id)
  if (length(matching)>0)
  {
    matching_purity <- purity[matching, data_list]
  }
  else
  {
    matching_purity <- rep(NA, length(data_list))
  }
  patient_purity <- rbind(patient_purity, matching_purity)
}
rownames(patient_purity) <- seq(1,nrow(patient_purity))
mode(patient_purity) <- 'numeric'
results_purity <- cbind(results, patient_purity)
results_purity$patient_purity <- as.numeric(as.character(results_purity$patient_purity))
results_all <- cbind(results_all, patient_purity)

write.table(results_all, "cdr3_results_with_expression_blood_ptprc_burden_purity_1078.txt", quote=FALSE, sep="\t", row.names=FALSE)

results_all <- results
results_all <- results_all[, !duplicated(colnames(results_all))]
cols_remove <- c("CD4", "CD8A", "FOXP3", "Q0VAE8", "P08575", "uc001gus.1", "P08575.2",
                 "F5GZM5", "uc001guv.1", "uc001guw.1", "Q5T9M4", "E9PKH0", "F5GXZ3")
results_all <- results_all[,!(names(results_all) %in% cols_remove)]
write.table(results_all, "cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_clustering_age_dedup_1078.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)
results <- results_all


### compare purity ###
results_purity <- results_purity[which(!is.na(results_purity$patient_purity)),]
split_subtypes <- split(results_purity, results_purity$clinical_group)
hr_positive <- split_subtypes[[1]]
her2_positive <- split_subtypes[[2]]
tnbc <- split_subtypes[[4]]

par(mfrow=c(2,2))
# purity
cdr3_rpm <- results_purity$exome_imseq/results_purity$exome_reads*1000000
purity <- results_purity$patient_purity
purity[which(purity>=median(purity))] <- 1
purity[which(purity<median(purity))] <- 0
purity_wilcox <- wilcox.test(cdr3_rpm ~ purity) # 2.2e-16
purity_wilcox
not <- cdr3_rpm[which(purity==0)]
pure <- cdr3_rpm[which(purity==1)]
boxplot(not, pure, names=c("Low", "High"), xlab="Purity", ylab="CDR3 RPM", main=c("BRCA", "n=1070"),
        ylim=c(0,0.14))
text(2, 0.12, "p=2.2e-16")
# purity by subtype
# her2+
cdr3_rpm <- her2_positive$exome_imseq/her2_positive$exome_reads*1000000
purity <- her2_positive$patient_purity
purity[which(purity>=median(purity))] <- 1
purity[which(purity<median(purity))] <- 0
purity_wilcox <- wilcox.test(cdr3_rpm ~ purity) # 2.2e-16
purity_wilcox
not <- cdr3_rpm[which(purity==0)]
pure <- cdr3_rpm[which(purity==1)]
boxplot(not, pure, names=c("Low", "High"), xlab="Purity", ylab="CDR3 RPM", main=c("HER2+", "n=157"),
        ylim=c(0,0.14))
text(2, 0.12, "p=0.002")
# hr+
cdr3_rpm <- hr_positive$exome_imseq/hr_positive$exome_reads*1000000
purity <- hr_positive$patient_purity
purity[which(purity>=median(purity))] <- 1
purity[which(purity<median(purity))] <- 0
purity_wilcox <- wilcox.test(cdr3_rpm ~ purity) # 2.2e-16
purity_wilcox
not <- cdr3_rpm[which(purity==0)]
pure <- cdr3_rpm[which(purity==1)]
boxplot(not, pure, names=c("Low", "High"), xlab="Purity", ylab="CDR3 RPM", main=c("HR+", "n=430"),
        ylim=c(0,0.14))
text(2, 0.12, "p=3.414e-09")
# tnbc
cdr3_rpm <- tnbc$exome_imseq/tnbc$exome_reads*1000000
purity <- tnbc$patient_purity
purity[which(purity>=median(purity))] <- 1
purity[which(purity<median(purity))] <- 0
purity_wilcox <- wilcox.test(cdr3_rpm ~ purity) # 2.2e-16
purity_wilcox
not <- cdr3_rpm[which(purity==0)]
pure <- cdr3_rpm[which(purity==1)]
boxplot(not, pure, names=c("Low", "High"), xlab="Purity", ylab="CDR3 RPM", main=c("TNBC", "n=113"),
        ylim=c(0,0.14))
text(2, 0.12, "p=0.017")
#
par(mfrow=c(1,1))

### compare clustering ###
results_clustering <- results_clustering[which(!is.na(results_clustering$CLUS_mRNAseq_cNMF)),]
split_subtypes <- split(results_clustering, results_clustering$clinical_group)
hr_positive <- split_subtypes[[1]]
her2_positive <- split_subtypes[[2]]
tnbc <- split_subtypes[[4]]

par(mfrow=c(2,2))
# clustering
cdr3_rpm <- results_clustering$exome_imseq/results_clustering$exome_reads*1000000
clustering <- as.character(results_clustering$CLUS_mRNAseq_cNMF)
cluster_1 <- cdr3_rpm[which(clustering==" 1")]
cluster_2 <- cdr3_rpm[which(clustering==" 2")]
cluster_3 <- cdr3_rpm[which(clustering==" 3")]
boxplot(cluster_1, cluster_2, cluster_3)#, names=c("Low", "High"), xlab="clustering", ylab="CDR3 RPM", main=c("BRCA", "n=1070"),
        #ylim=c(0,0.14))
text(2, 0.12, "p=2.2e-16")
# clustering by subtype
# her2+
cdr3_rpm <- her2_positive$exome_imseq/her2_positive$exome_reads*1000000
clustering <- her2_positive$patient_clustering
clustering[which(clustering>=median(clustering))] <- 1
clustering[which(clustering<median(clustering))] <- 0
clustering_wilcox <- wilcox.test(cdr3_rpm ~ clustering) # 2.2e-16
clustering_wilcox
not <- cdr3_rpm[which(clustering<median(clustering))]
pure <- cdr3_rpm[which(clustering>=median(clustering))]
boxplot(not, pure, names=c("Low", "High"), xlab="clustering", ylab="CDR3 RPM", main=c("HER2+", "n=157"),
        ylim=c(0,0.14))
text(2, 0.12, "p=0.002")
# hr+
cdr3_rpm <- hr_positive$exome_imseq/hr_positive$exome_reads*1000000
clustering <- hr_positive$patient_clustering
clustering[which(clustering>=median(clustering))] <- 1
clustering[which(clustering<median(clustering))] <- 0
clustering_wilcox <- wilcox.test(cdr3_rpm ~ clustering) # 2.2e-16
clustering_wilcox
not <- cdr3_rpm[which(clustering<median(clustering))]
pure <- cdr3_rpm[which(clustering>=median(clustering))]
boxplot(not, pure, names=c("Low", "High"), xlab="clustering", ylab="CDR3 RPM", main=c("HR+", "n=430"),
        ylim=c(0,0.14))
text(2, 0.12, "p=3.414e-09")
# tnbc
cdr3_rpm <- tnbc$exome_imseq/tnbc$exome_reads*1000000
clustering <- tnbc$patient_clustering
clustering[which(clustering>=median(clustering))] <- 1
clustering[which(clustering<median(clustering))] <- 0
clustering_wilcox <- wilcox.test(cdr3_rpm ~ clustering) # 2.2e-16
clustering_wilcox
not <- cdr3_rpm[which(clustering<median(clustering))]
pure <- cdr3_rpm[which(clustering>=median(clustering))]
boxplot(not, pure, names=c("Low", "High"), xlab="clustering", ylab="CDR3 RPM", main=c("TNBC", "n=113"),
        ylim=c(0,0.14))
text(2, 0.12, "p=0.017")
#
par(mfrow=c(1,1))

### compare mutation ###
results_mutation <- results_mutation[which(!is.na(results_mutation$rate_non)),]
split_subtypes <- split(results_mutation, results_mutation$clinical_group)
hr_positive <- split_subtypes[[1]]
her2_positive <- split_subtypes[[2]]
tnbc <- split_subtypes[[4]]

par(mfrow=c(2,2))
# mutation
cdr3_rpm <- results_mutation$exome_imseq/results_mutation$exome_reads*1000000
#mutation <- as.numeric(as.character(results_mutation$rate_non))
mutation <- as.numeric(as.character(results_mutation$rate_non)) + as.numeric(as.character(results_mutation$rate_sil))
mutation[which(mutation>=median(mutation))] <- 1
mutation[which(mutation<median(mutation))] <- 0
mutation_wilcox <- wilcox.test(cdr3_rpm ~ mutation) # 2.2e-16
mutation_wilcox
not <- cdr3_rpm[which(mutation==0)]
pure <- cdr3_rpm[which(mutation==1)]
boxplot(not, pure, names=c("Low", "High"), xlab="mutation", ylab="CDR3 RPM", main=c("BRCA", "n=1070"),
        ylim=c(0,0.14))
text(2, 0.12, "p=2.2e-16")
# mutation by subtype
# her2+
cdr3_rpm <- her2_positive$exome_imseq/her2_positive$exome_reads*1000000
#mutation <- as.numeric(as.character(her2_positive$rate_non))
mutation <- as.numeric(as.character(her2_positive$rate_non)) + as.numeric(as.character(her2_positive$rate_sil))
mutation[which(mutation>=median(mutation))] <- 1
mutation[which(mutation<median(mutation))] <- 0
mutation_wilcox <- wilcox.test(cdr3_rpm ~ mutation) # 2.2e-16
mutation_wilcox
not <- cdr3_rpm[which(mutation<median(mutation))]
pure <- cdr3_rpm[which(mutation>=median(mutation))]
boxplot(not, pure, names=c("Low", "High"), xlab="mutation", ylab="CDR3 RPM", main=c("HER2+", "n=157"),
        ylim=c(0,0.14))
text(2, 0.12, "p=0.002")
# hr+
cdr3_rpm <- hr_positive$exome_imseq/hr_positive$exome_reads*1000000
#mutation <- as.numeric(as.character(results_mutation$rate_non))
mutation <- as.numeric(as.character(hr_positive$rate_non)) + as.numeric(as.character(hr_positive$rate_sil))
mutation[which(mutation>=median(mutation))] <- 1
mutation[which(mutation<median(mutation))] <- 0
mutation_wilcox <- wilcox.test(cdr3_rpm ~ mutation) # 2.2e-16
mutation_wilcox
not <- cdr3_rpm[which(mutation<median(mutation))]
pure <- cdr3_rpm[which(mutation>=median(mutation))]
boxplot(not, pure, names=c("Low", "High"), xlab="mutation", ylab="CDR3 RPM", main=c("HR+", "n=430"),
        ylim=c(0,0.14))
text(2, 0.12, "p=3.414e-09")
# tnbc
cdr3_rpm <- tnbc$exome_imseq/tnbc$exome_reads*1000000
#mutation <- as.numeric(as.character(results_mutation$rate_non))
mutation <- as.numeric(as.character(tnbc$rate_non)) + as.numeric(as.character(tnbc$rate_sil))
mutation[which(mutation>=median(mutation))] <- 1
mutation[which(mutation<median(mutation))] <- 0
mutation_wilcox <- wilcox.test(cdr3_rpm ~ mutation) # 2.2e-16
mutation_wilcox
not <- cdr3_rpm[which(mutation<median(mutation))]
pure <- cdr3_rpm[which(mutation>=median(mutation))]
boxplot(not, pure, names=c("Low", "High"), xlab="mutation", ylab="CDR3 RPM", main=c("TNBC", "n=113"),
        ylim=c(0,0.14))
text(2, 0.12, "p=0.017")
#
par(mfrow=c(1,1))