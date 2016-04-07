library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)

#results <- read.csv("cdr3_results_with_groups_complete.txt", sep="\t")
#results <- read.csv("cdr3_results_with_expression.txt", sep="\t")
results <- read.csv("/Users/Eric/cdr3_results_with_blood_1078.txt", sep="\t")
idxstats <- "/Users/Eric/tcga/clonotypes/final/rna_idxstats/"

exome_rpm <- as.numeric(results$exome_imseq)/as.numeric(results$exome_reads)*1000000
rna_rpm <- as.numeric(results$rna_imseq)/as.numeric(results$rna_reads)*1000000

# get idxstats values
file_summaries_rna <- list.files(path=idxstats, pattern="*.txt", full.names=T, recursive=FALSE)
tcrb <- c()
for (i in 1:nrow(results))
{
  if (is.na(results$rna_id[i])) {tcrb <- c(tcrb, NA)}
  else
  {
    idx_file <- paste(results$rna_id[i], ".TCR.sorted.bam.idxstats.txt", sep="")
    idx_path <- paste(idxstats, idx_file, sep="")
    if (file.exists(idx_path))
    {
      idx <- read.csv(idx_path, sep="\t", header=FALSE, row.names=1)
      tcrb <- c(tcrb, idx["chr7",2])         
    }
    else {tcrb <- c(tcrb, NA)}
  }

}

tcrb_rpm <- as.numeric(tcrb)/as.numeric(results$rna_reads)*1000000
results <- cbind(results, tcrb)

# change cdr3 cutoff to nonzero
cdr3_nonzero <- c()
for (i in 1:nrow(results))
{
  if (results$exome_imseq[i]>0) {cdr3_nonzero <- c(cdr3_nonzero, "high")}
  else {cdr3_nonzero <- c(cdr3_nonzero, "low")}
}
results <- cbind(results, cdr3_nonzero)

# change old cdr3 cutoff to two groups
#cdr3_two <- c()
#for (i in 1:nrow(results))
#{
#  if (as.character(results$cdr3_grouping[i])=="zero") {cdr3_two <- c(cdr3_two, "low")}
#  else {cdr3_two <- c(cdr3_two, as.character(results$cdr3_grouping[i]))}
#}
#results <- cbind(results, cdr3_two)

# redo lymphocyte cutoff global
lymphocyte_presence <- c()
lymph_med <- median(results$lymphocyte_percent)
for (i in 1:nrow(results))
{
  lymph <- results$lymphocyte_percent[i]
  if (lymph <= lymph_med) {lymphocyte_presence <- c(lymphocyte_presence, "0")}
  else {lymphocyte_presence <- c(lymphocyte_presence, "1")}
}
results$lymphocyte_presence <- as.numeric(lymphocyte_presence)

# make some comparisons
split_subtypes <- split(results, results$clinical_group)
hr_positive <- split_subtypes[[1]]
her2_positive <- split_subtypes[[2]]
tnbc <- split_subtypes[[4]]

# redo lymphocyte cutoff by subtype
lymphocyte_presence <- c()
lymph_med <- median(hr_positive$lymphocyte_percent)
for (i in 1:nrow(hr_positive))
{
  lymph <- hr_positive$lymphocyte_percent[i]
  if (lymph <= lymph_med) {lymphocyte_presence <- c(lymphocyte_presence, "0")}
  else {lymphocyte_presence <- c(lymphocyte_presence, "1")}
}
hr_positive$lymphocyte_presence <- as.numeric(lymphocyte_presence)

lymphocyte_presence <- c()
lymph_med <- median(her2_positive$lymphocyte_percent)
for (i in 1:nrow(her2_positive))
{
  lymph <- her2_positive$lymphocyte_percent[i]
  if (lymph <= lymph_med) {lymphocyte_presence <- c(lymphocyte_presence, "0")}
  else {lymphocyte_presence <- c(lymphocyte_presence, "1")}
}
her2_positive$lymphocyte_presence <- as.numeric(lymphocyte_presence)

lymphocyte_presence <- c()
lymph_med <- median(tnbc$lymphocyte_percent)
for (i in 1:nrow(tnbc))
{
  lymph <- tnbc$lymphocyte_percent[i]
  if (lymph <= lymph_med) {lymphocyte_presence <- c(lymphocyte_presence, "0")}
  else {lymphocyte_presence <- c(lymphocyte_presence, "1")}
}
tnbc$lymphocyte_presence <- as.numeric(lymphocyte_presence)

# survival by cdr3
pdf("hr_positive_cdr3_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv <- hr_positive[,c("cdr3_nonzero", "days", "vital")]
surv$days <- surv$days/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", xlim=c(0,11),
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2) 
#legend(x="bottomleft", c("nonzero", "zero"), fill=c(1,2), cex=2, bty="n")
legend(x="bottomleft", c("iDNA>0", "iDNA=0"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1]) # p=0.412, n=433, events=20
text(x=9.25, y=0.95, labels="p=0.415", cex=2) # changed text to HR p-value
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

pdf("her2_positive_cdr3_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv <- her2_positive[,c("cdr3_nonzero", "days", "vital")]
surv$days <- surv$days/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", xlim=c(0,11),
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2) 
#legend(x="bottomleft", c("nonzero", "zero"), fill=c(1,2), cex=2, bty="n")
legend(x="bottomleft", c("iDNA>0", "iDNA=0"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1]) # p=0.016, n=157, events=18
text(x=9.5, y=0.95, labels="p=0.022", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

pdf("tnbc_cdr3_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv <- tnbc[,c("cdr3_nonzero", "days", "vital")]
surv$days <- surv$days/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", xlim=c(0,11),
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2) 
#legend(x="bottomleft", c("nonzero", "zero"), fill=c(1,2), cex=2, bty="n")
legend(x="bottomleft", c("iDNA>0", "iDNA=0"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1]) # p=0.057, n=114, events=13
text(x=9.5, y=0.95, labels="p=0.071", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

# survival by lymph
pdf("hr_positive_til_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv <- hr_positive[,c("lymphocyte_presence", "days", "vital")]
surv$days <- surv$days/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(2,1), xlab="Survival Time (years)", ylab="Survival Probability", xlim=c(0,11),
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2) 
legend(x="bottomleft", c("TIL high", "TIL low"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1]) # p=0.947
text(x=9.5, y=0.95, labels="p=0.947", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

pdf("her2_positive_til_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv <- her2_positive[,c("lymphocyte_presence", "days", "vital")]
surv$days <- surv$days/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(2,1), xlab="Survival Time (years)", ylab="Survival Probability", xlim=c(0,11),
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2) 
legend(x="bottomleft", c("TIL high", "TIL low"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1]) # p=0.166
text(x=9.5, y=0.95, labels="p=0.176", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

pdf("tnbc_til_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv <- tnbc[,c("lymphocyte_presence", "days", "vital")]
surv$days <- surv$days/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(2,1), xlab="Survival Time (years)", ylab="Survival Probability", xlim=c(0,11),
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2) 
legend(x="bottomleft", c("TIL high", "TIL low"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1]) # p=0.805
text(x=9.5, y=0.95, labels="p=0.805", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

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
genes <- c("CD4", "CD8A", "FOXP3")
cd4 <- c()
cd8a <- c()
foxp3 <- c()
for (i in 1:nrow(results))
{
  patient <- as.character(results$patient_id[i])
  matching <- which(rownames(expression)==patient)
  if (length(matching)>0)
  {
    matching_expression <- expression[matching, genes]
    cd4 <- c(cd4, expression[matching, "CD4"])
    cd8a <- c(cd8a, expression[matching, "CD8A"])
    foxp3 <- c(foxp3, expression[matching, "FOXP3"])
  }
  else
  {
    cd4 <- c(cd4, NA)
    cd8a <- c(cd8a, NA)
    foxp3 <- c(foxp3, NA)
  }
}
results <- cbind(results, cd4, cd8a, foxp3)

# make some comparisons
split_subtypes <- split(results, results$clinical_group)
hr_positive <- split_subtypes[[1]]
her2_positive <- split_subtypes[[2]]
tnbc <- split_subtypes[[4]]

# cd4 exome
plot(hr_positive$exome_rpm, hr_positive$cd4)
plot(her2_positive$exome_rpm, her2_positive$cd4)
plot(tnbc$exome_rpm, tnbc$cd4)

cor.test(hr_positive$exome_rpm, hr_positive$cd4)
cor.test(her2_positive$exome_rpm, her2_positive$cd4)
cor.test(tnbc$exome_rpm, tnbc$cd4)

# rna
plot(hr_positive$rna_rpm, hr_positive$cd4)
plot(her2_positive$rna_rpm, her2_positive$cd4)
plot(tnbc$rna_rpm, tnbc$cd4)

cor.test(hr_positive$rna_rpm, hr_positive$cd4)
cor.test(her2_positive$rna_rpm, her2_positive$cd4)
cor.test(tnbc$rna_rpm, tnbc$cd4)

# tcrb
plot(hr_positive$rna_rpm, hr_positive$tcrb_rpm)
plot(her2_positive$rna_rpm, her2_positive$tcrb_rpm)
plot(tnbc$rna_rpm, tnbc$tcrb_rpm)

cor.test(hr_positive$rna_rpm, hr_positive$tcrb_rpm)
cor.test(her2_positive$rna_rpm, her2_positive$tcrb_rpm)
cor.test(tnbc$rna_rpm, tnbc$tcrb_rpm)

### correlation table output ###
all_correlations <- function(cdr3_table)
{
  output_table <- matrix(nrow=6, ncol=4)
  rownames(output_table) <- c("Exome_P", "RNA_P", "TIL_P", "Exome_S", "RNA_S", "TIL_S")
  colnames(output_table) <- c("CD4", "CD8A", "FOXP3", "TCRB")
  output_table[1,1] <- cor.test(cdr3_table$exome_rpm, cdr3_table$cd4, method="pearson")[[3]]
  output_table[1,2] <- cor.test(cdr3_table$exome_rpm, cdr3_table$cd8a, method="pearson")[[3]]
  output_table[1,3] <- cor.test(cdr3_table$exome_rpm, cdr3_table$foxp3, method="pearson")[[3]]
  output_table[1,4] <- cor.test(cdr3_table$exome_rpm, cdr3_table$tcrb_rpm, method="pearson")[[3]]
  output_table[2,1] <- cor.test(cdr3_table$rna_rpm, cdr3_table$cd4, method="pearson")[[3]]
  output_table[2,2] <- cor.test(cdr3_table$rna_rpm, cdr3_table$cd8a, method="pearson")[[3]]
  output_table[2,3] <- cor.test(cdr3_table$rna_rpm, cdr3_table$foxp3, method="pearson")[[3]]
  output_table[2,4] <- cor.test(cdr3_table$rna_rpm, cdr3_table$tcrb_rpm, method="pearson")[[3]]
  output_table[3,1] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$cd4, method="pearson")[[3]]
  output_table[3,2] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$cd8a, method="pearson")[[3]]
  output_table[3,3] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$foxp3, method="pearson")[[3]]
  output_table[3,4] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$tcrb_rpm, method="pearson")[[3]]
  output_table[4,1] <- cor.test(cdr3_table$exome_rpm, cdr3_table$cd4, method="spearman")[[3]]
  output_table[4,2] <- cor.test(cdr3_table$exome_rpm, cdr3_table$cd8a, method="spearman")[[3]]
  output_table[4,3] <- cor.test(cdr3_table$exome_rpm, cdr3_table$foxp3, method="spearman")[[3]]
  output_table[4,4] <- cor.test(cdr3_table$exome_rpm, cdr3_table$tcrb_rpm, method="spearman")[[3]]
  output_table[5,1] <- cor.test(cdr3_table$rna_rpm, cdr3_table$cd4, method="spearman")[[3]]
  output_table[5,2] <- cor.test(cdr3_table$rna_rpm, cdr3_table$cd8a, method="spearman")[[3]]
  output_table[5,3] <- cor.test(cdr3_table$rna_rpm, cdr3_table$foxp3, method="spearman")[[3]]
  output_table[5,4] <- cor.test(cdr3_table$rna_rpm, cdr3_table$tcrb_rpm, method="spearman")[[3]]
  output_table[6,1] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$cd4, method="spearman")[[3]]
  output_table[6,2] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$cd8a, method="spearman")[[3]]
  output_table[6,3] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$foxp3, method="spearman")[[3]]
  output_table[6,4] <- cor.test(cdr3_table$lymphocyte_percent, cdr3_table$tcrb_rpm, method="spearman")[[3]]
  return(output_table)
}

hr_positive_cors <- all_correlations(hr_positive)
her2_positive_cors <- all_correlations(her2_positive)
tnbc_cors <- all_correlations(tnbc)

#results_out <- results[,c(1:18,42,45:47)]
#colnames(results_out)[19] <- "tcrb_reads"
#write.table(results_out, "cdr3_results_with_expression.txt", quote=FALSE, sep="\t", row.names=FALSE)

colnames(results)[24] <- "tcrb_reads"
results <- results[,-25]
write.table(results, "cdr3_results_with_expression_blood_1078.txt", quote=FALSE, sep="\t", row.names=FALSE)
