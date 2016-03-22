library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)

### messing around with lymphocytes
results_all <- read.csv("cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_1078.txt")
results <- results_all[!is.na(results_all$T.regs),]

# get rid of her2+
#results <- results[which(results$clinical_group!="her2+"),]
# get rid of "other"
#results <- results[which(results$clinical_group!="other"),]

# isolate gene signature values
#signatures <- results[,20:41] # old cdr3_results_with_gsva
signatures <- results[,50:71] # new cdr3_results_with lots more
rownames(signatures) <- results[,1]

# normalize?
signatures <- scale(signatures)

# Determine number of clusters
wss <- (nrow(signatures)-1)*sum(apply(signatures,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(signatures, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(signatures, 4)
# get cluster means 
aggregate(signatures,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(results, fit$cluster)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
clusplot(signatures, fit$cluster, color=TRUE, shade=TRUE, labels=0, lines=0)
# Centroid Plot against 1st 2 discriminant functions
plotcluster(signatures, fit$cluster)

for (i in 1:nrow(signatures))
{
  if (i%%100!=0) {rownames(signatures)[i] <- i}
}

# hierarchical cluster analysis
d <- dist(as.matrix(signatures))
hc <- hclust(d)
plot(hc, labels=FALSE)
mycl <- cutree(hc, h=9)
#mycl <- cutree(hc, h=max(hc$height/1.5))

results_hc <- cbind(results, mycl)

# add clusters to the 1078 results table (6 without the expression data)
mycl_1078 <- c()
names(mycl) <- results_hc[,1]
for (i in 1:nrow(results_all))
{
  patient_id <- as.character(results_all$patient_id[i])
  mycl_1078 <- c(mycl_1078, mycl[patient_id])
}
rownames(mycl_1078) <- results_all$patient_id
results_hc_1078 <- cbind(results_all, mycl_1078)
colnames(results_hc_1078)[72] <- "gsva_cluster"
write.table(results_hc_1078, "cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_clustering_four_groups_1078.txt", sep="\t", quote=FALSE, row.names=FALSE)

col_mat <- matrix(nrow=nrow(results), ncol=3)
exome_freq <- results$exome_imseq/results$exome_reads
rna_freq <- results$rna_imseq/results$rna_reads
exome_med <- median(exome_freq[which(exome_freq!=0)])
rna_med <- median(rna_freq[which(rna_freq!=0)])
for (i in 1:nrow(col_mat))
{
  group <- results$clinical_group[i]
  exome <- exome_freq[i]
  rna <- rna_freq[i]
  if (exome<=exome_med) {col_mat[i,2]="white"}
  else {col_mat[i,2]="black"}
  if (rna<=rna_med) {col_mat[i,3]="white"}
  else {col_mat[i,3]="purple"}
  if (group=="her2-") {col_mat[i,1]="red"}
  else {col_mat[i,1]="green"}
}

heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none")

# survival of tnbc clusters
results <- read.csv("cdr3_results_with_gsva.csv")
results <- na.omit(results)
results <- results[which(results$clinical_group=="triple-"),]
# cluster gene signatures
signatures <- results[,20:41]
rownames(signatures) <- results[,1]
signatures <- scale(signatures)
d <- dist(as.matrix(signatures))
hc <- hclust(d)
mycl <- cutree(hc, h=10)

vital <- results$vital
d <- results$days/365
surv <- cbind(mycl, d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])

plot(my.fit, col=c(1,2), main="Survival of tnbc by signature clusters", xlab="Survival Time (years)", ylab="Survival Probability")
#legend(0.5,0.4, c("high", "low", "zero"), fill=c(1,2,3))

# survival of her2+ clusters
results <- read.csv("cdr3_results_with_gsva.csv")
results <- na.omit(results)
results <- results[which(results$clinical_group=="her2+"),]
# cluster gene signatures
signatures <- results[,20:41]
rownames(signatures) <- results[,1]
signatures <- scale(signatures)
d <- dist(as.matrix(signatures))
hc <- hclust(d)
mycl <- cutree(hc, h=11)

vital <- results$vital
d <- results$days/365
surv <- cbind(mycl, d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])

plot(my.fit, col=c(1,2), main="Survival of her2+ by signature clusters", xlab="Survival Time (years)", ylab="Survival Probability")
#legend(0.5,0.4, c("high", "low", "zero"), fill=c(1,2,3))

survdiff(my.surv ~ surv[,1])

# certain metagenes by subtype
# TNBC
results <- read.csv("cdr3_results_with_gsva.csv")
results <- na.omit(results)
tnbc <- results[which(results$clinical_group=="triple-"),]
masta <- tnbc$Mast.a
treg <- tnbc$T.regs
macro2 <- tnbc$Macro.2

signature_levels <- c()
for (i in 1:nrow(tnbc))
{
  if (masta[i] > median(masta)) {masta_row <- "high"}
  else {masta_row <- "low"}
  if (treg[i] > median(treg)) {treg_row <- "high"}
  else {treg_row <- "low"}
  if (macro2[i] > median(macro2)) {macro2_row <- "high"}
  else {macro2_row <- "low"}
  signature_levels <- rbind(signature_levels, c(masta_row, treg_row, macro2_row))
}
colnames(signature_levels) <- c("Mast.a_group", "T.Reg_group", "Macro.2_group")

vital <- tnbc$vital
d <- tnbc$days/365
surv <- cbind(signature_levels[,1], d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of tnbc by masta", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low"), fill=c(1,2))


vital <- tnbc$vital
d <- tnbc$days/365
surv <- cbind(signature_levels[,2], d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of tnbc by treg", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low"), fill=c(1,2))

vital <- tnbc$vital
d <- tnbc$days/365
surv <- cbind(signature_levels[,3], d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of tnbc by macro2", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low"), fill=c(1,2))

# her2+
results <- read.csv("cdr3_results_with_gsva.csv")
results <- na.omit(results)
her2 <- results[which(results$clinical_group=="her2+"),]
masta <- her2$Mast.a
treg <- her2$T.regs
macro2 <- her2$Macro.2

signature_levels <- c()
for (i in 1:nrow(her2))
{
  if (masta[i] > median(masta)) {masta_row <- "high"}
  else {masta_row <- "low"}
  if (treg[i] > median(treg)) {treg_row <- "high"}
  else {treg_row <- "low"}
  if (macro2[i] > median(macro2)) {macro2_row <- "high"}
  else {macro2_row <- "low"}
  signature_levels <- rbind(signature_levels, c(masta_row, treg_row, macro2_row))
}
colnames(signature_levels) <- c("Mast.a_group", "T.Reg_group", "Macro.2_group")

vital <- her2$vital
d <- her2$days/365
surv <- cbind(signature_levels[,1], d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of her2+ by masta", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low"), fill=c(1,2))
survdiff(my.surv ~ surv[,1])

vital <- her2$vital
d <- her2$days/365
surv <- cbind(signature_levels[,2], d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of her2+ by treg", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low"), fill=c(1,2))

vital <- her2$vital
d <- her2$days/365
surv <- cbind(signature_levels[,3], d, vital)
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of her2+ by macro2", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low"), fill=c(1,2))

# stacked boxplot by subtype
results <- read.csv("cdr3_results_with_gsva.csv")
results <- na.omit(results)
results <- results[which(results$clinical_group!="other"),]
results$clinical_group <- factor(results$clinical_group)
counts <- table(results$cdr3_grouping, results$clinical_group)
counts_prop <- prop.table(counts, 2)
line1 <- colnames(counts_prop)
line1[1] <- "hr+"
line2 <- colSums(counts)
names_arg <- paste(line1, line2, sep = "\n")
colnames(counts_prop) <- names_arg
barplot(counts_prop, main="Distribution of CDR3 groups", xlab="BRCA Subtype",
        ylab="Percentage of Patients", col=c("blue", "red", "yellow"),
        legend=rownames(counts))

# survival plots by cdr3 low/high within subtype
results <- read.csv("cdr3_results_with_gsva.csv")
results_split <- split(results, results$clinical_group)
hr <- results_split[[1]]
her2 <- results_split[[2]]
tnbc <- results_split[[4]]

replace_cdr3_grouping <- function(subtype)
{
  subtype_freq <- subtype$exome_imseq/subtype$exome_reads
  subtype_med <- median(subtype_freq[which(subtype_freq!=0)])
  for (i in 1:nrow(subtype))
  {
    row_freq <- subtype_freq[i]
    if (row_freq == 0) {subtype$cdr3_grouping[i] <- "zero"}
    else {subtype$cdr3_grouping[i] <- "high"}
  }
  return(subtype)
}

results_hc <- replace_cdr3_grouping(hr)
results_her2 <- replace_cdr3_grouping(her2)
results_tnbc <- replace_cdr3_grouping(tnbc)

surv_hr <- results_hc[,c("cdr3_grouping", "days", "vital")]
surv_hr$days <- results_hc$days/365
my.surv_hr <- Surv(as.numeric(surv_hr[,2]), as.numeric(surv_hr[,3]))
my.fit_hr <- survfit(my.surv_hr ~ surv_hr[,1])
plot(my.fit_hr, col=c(1,2), main="Survival of HR+ by CDR3 groups", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("nonzero", "zero"), fill=c(1,2), cex=0.75)

surv_her2 <- results_her2[,c("cdr3_grouping", "days", "vital")]
surv_her2$days <- results_her2$days/365
my.surv_her2 <- Surv(as.numeric(surv_her2[,2]), as.numeric(surv_her2[,3]))
my.fit_her2 <- survfit(my.surv_her2 ~ surv_her2[,1])
plot(my.fit_her2, col=c(1,2), main="Survival of HER2+ by CDR3 groups", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("nonzero", "zero"), fill=c(1,2), cex=0.75)

surv_tnbc <- results_tnbc[,c("cdr3_grouping", "days", "vital")]
surv_tnbc$days <- results_tnbc$days/365
my.surv_tnbc <- Surv(as.numeric(surv_tnbc[,2]), as.numeric(surv_tnbc[,3]))
my.fit_tnbc <- survfit(my.surv_tnbc ~ surv_tnbc[,1])
plot(my.fit_tnbc, col=c(1,2), main="Survival of TNBC+ by CDR3 groups", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("nonzero", "zero"), fill=c(1,2), cex=0.75)

#survdiff(my.surv_clusters ~ surv_clusters[,1])