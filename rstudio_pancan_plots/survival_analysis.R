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

exome_rpm_med <- median(filter(cdr3_results_current, exome_level>0)$exome_rpm)
cdr3_results_current<-mutate(cdr3_results_current,exome_lowhigh=ifelse(exome_rpm>exome_rpm_med,"1","0"))

til_med <- median(cdr3_results_current$percent_tils, na.rm=TRUE)
cdr3_results_current<-mutate(cdr3_results_current,til_lowhigh=ifelse(percent_tils>til_med,"1","0"))

cdr3_results_current<-mutate(cdr3_results_current,gsva_lowhigh=ifelse(gsva_cluster==1 | gsva_cluster==4,"1","0"))

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

# survival
surv_results_p <- function(surv_matrix)
{
  p_val_out <- rep(NA, 5)
  # exome nonzero
  surv <- surv_matrix[,c("exome_level", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  if (nrow(my.fit)>1) {p_val_out[1] <- summary(coxph(my.surv ~ surv[,1], method="breslow"))$logtest[3]}

  # RNA low/high
  surv <- surv_matrix[,c("RNA_level", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  if (nrow(my.fit)>1) {p_val_out[2] <- summary(coxph(my.surv ~ surv[,1], method="breslow"))$logtest[3]}
  
  # exome low/high
  surv <- surv_matrix[,c("exome_lowhigh", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  if (nrow(my.fit)>1) {p_val_out[3] <- summary(coxph(my.surv ~ surv[,1], method="breslow"))$logtest[3]}
  
  # TIL low/high
  surv <- surv_matrix[,c("til_lowhigh", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  if (nrow(my.fit)>1) {p_val_out[4] <- summary(coxph(my.surv ~ surv[,1], method="breslow"))$logtest[3]}
  
  # GSVA low/high
  surv <- surv_matrix[,c("gsva_lowhigh", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
  if (nrow(my.fit)>1) {p_val_out[5] <- summary(coxph(my.surv ~ surv[,1], method="breslow"))$logtest[3]}
  
  return(p_val_out)
}

types <- split(cancer_current, cancer_current$cohort)
#types_surv_exome <- lapply(types, surv_results_exome)
#types_surv_RNA <- lapply(types, surv_results_RNA)
types_surv_p <- sapply(types, surv_results_p)
rownames(types_surv_p) <- c("exome zero/nonzero", "rna low/high", "exome low/high",
                            "TILs low/high", "GSVA low/high")

write.table(types_surv_p, "survival_pvals.txt", quote=FALSE, sep="\t")

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

### FIGURE: BRCA survival by GSVA low/high ###
pdf("brca_gsva_lowhigh_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$BRCA
surv <- surv_matrix[,c("gsva_lowhigh", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="bottomright", c("Low", "High"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.04", cex=2)
text(x=5, y=0.2, labels="n=72/32 low/high", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

### FIGURE: HNSC survival by GSVA low/high ###
pdf("hnsc_gsva_lowhigh_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$HNSC
surv <- surv_matrix[,c("gsva_lowhigh", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="topright", c("Low", "High"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.04", cex=2)
text(x=5, y=0.1, labels="n=103/64 low/high", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)


### FIGURE: CESC survival by GSVA low/high ###
pdf("cesc_gsva_lowhigh_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$CESC
surv <- surv_matrix[,c("gsva_lowhigh", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="topright", c("Low", "High"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.02", cex=2)
text(x=5, y=0.1, labels="n=40/20 low/high", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)


### FIGURE: LGG survival by GSVA low/high ###
pdf("lgg_gsva_lowhigh_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$LGG
surv <- surv_matrix[,c("gsva_lowhigh", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="topright", c("Low", "High"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.01", cex=2)
text(x=5, y=0.1, labels="n=58/34 low/high", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)

### FIGURE: LIHC survival by GSVA low/high ###
pdf("lihc_gsva_lowhigh_survival.pdf", width=8, height=8)
mar.default <- c(5.1, 4.1, 4.1, 2.1)
oma.default <- c(0,0,0,0)
par(mar = mar.default + c(1,1,0,0))
par(oma = oma.default + c(0,1,0,0))
surv_matrix <- types$LIHC
surv <- surv_matrix[,c("gsva_lowhigh", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), xlab="Survival Time (years)", ylab="Survival Probability", 
     cex.axis=1.75, cex.lab=2.25, cex.main=1.75, lwd=4, cex=2)
legend(x="topright", c("Low", "High"), fill=c(1,2), cex=2, bty="n")
survdiff(my.surv ~ surv[,1])
text(x=9.5, y=0.95, labels="p=0.01", cex=2)
text(x=5, y=0.1, labels="n=65/24 low/high", cex=2)
dev.off()
cox.fit <- coxph(my.surv ~ surv[,1])
summary(cox.fit)