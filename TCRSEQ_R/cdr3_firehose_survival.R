library(OIsurv)
library(cgdsr)

### messing around with lymphocytes
results <- read.csv("cdr3_results_with_gsva.csv")
#results <- results[,-c(15:17)]  # remove the er, pr, her2
results <- results[complete.cases(results),]
exome_freq <- as.numeric(results$exome_imseq)/as.numeric(results$exome_reads)
rna_freq <- as.numeric(results$rna_imseq)/as.numeric(results$rna_reads)
treg_sig <- results$T.regs
lymphocytes <- results$lymphocyte_percent
exome_cdr3_2 <- rep("",nrow(results))
rna_cdr3_2 <- rep("",nrow(results))
exome_cdr3_3 <- rep("",nrow(results))
rna_cdr3_3 <- rep("",nrow(results))
lymph_group <- rep("",nrow(results))
exome_rna <- rep("",nrow(results))
treg_group <- rep("",nrow(results))
exome_treg <- rep("",nrow(results))
for (i in 1:nrow(results))
{
  exome <- exome_freq[i]
  rna <- rna_freq[i]
  lymph <- lymphocytes[i]
  treg <- treg_sig[i]
  if (exome==0)
  {
    exome_cdr3_2[i] <- "low"
    exome_cdr3_3[i] <- "zero"
  }
  else if (exome<=median(exome_freq[which(exome_freq!=0)])) 
  {
    exome_cdr3_2[i] <- "low"
    exome_cdr3_3[i] <- "low"
  }
  else 
  {
    exome_cdr3_2[i] <- "high"
    exome_cdr3_3[i] <- "high"
  }
  
  if (rna==0)
  {
    rna_cdr3_2[i] <- "low"
    rna_cdr3_3[i] <- "zero"
  }
  else if (rna<=median(rna_freq[which(rna_freq!=0)])) 
  {
    rna_cdr3_2[i] <- "low"
    rna_cdr3_3[i] <- "low"
  }
  else 
  {
    rna_cdr3_2[i] <- "high"
    rna_cdr3_3[i] <- "high"
  }
  
  if (lymph<median(lymphocytes)) {lymph_group[i] <- "low"}
  else {lymph_group[i] <- "high"}
  
  if (exome==0 && rna==0) {exome_rna[i] <- "D-R-"}
  else if (exome==0 && rna>0) {exome_rna[i] <- "D-R+"}
  else if (exome>0 && rna==0) {exome_rna[i] <- "D+R-"}
  else {exome_rna[i] <- "D+R+"}
  
  if (treg<=median(treg_sig)) {treg_group[i] <- "low"}
  else {treg_group[i] <- "high"}
  
  if (exome_cdr3_2[i]=="high" && treg_group[i]=="high") {exome_treg[i] <- "T-high/CDR3-high"}
  else if (exome_cdr3_2[i]=="low" && treg_group[i]=="high") {exome_treg[i] <- "T-high/CDR3-low"}
  else if (exome_cdr3_2[i]=="high" && treg_group[i]=="low") {exome_treg[i] <- "T-low/CDR3-high"}
  else {exome_treg[i] <- "T-low/CDR3-low"}
}

vital <- results$vital
d <- results$days/365

# test
#survdiff(my.surv1 ~ surv_complete_lymph[,1])
#survdiff(my.surv2 ~ surv_complete_cdr3[,1])

groupings <- cbind(clinical_group, exome_cdr3_2, exome_cdr3_3, rna_cdr3_2, rna_cdr3_3, exome_rna, exome_treg, d, vital)
groupings <- as.data.frame(groupings)

# 3 group exome cdr3 by subtype
split_subtypes <- split(groupings, clinical_group)
hr_positive <- split_subtypes[[1]]
her2_positive <- split_subtypes[[2]]
tnbc <- split_subtypes[[4]]

surv_her2_positive_exome_3 <- as.matrix(her2_positive[,c("exome_cdr3_3", "d", "vital")])
my.surv_her2_positive_exome_3 <- Surv(as.numeric(surv_her2_positive_exome_3[,2]), as.numeric(surv_her2_positive_exome_3[,3]))
my.fit_her2_positive_exome_3 <- survfit(my.surv_her2_positive_exome_3 ~ surv_her2_positive_exome_3[,1])

surv_hr_positive_exome_3 <- as.matrix(hr_positive[,c("exome_cdr3_3", "d", "vital")])
my.surv_hr_positive_exome_3 <- Surv(as.numeric(surv_hr_positive_exome_3[,2]), as.numeric(surv_hr_positive_exome_3[,3]))
my.fit_hr_positive_exome_3 <- survfit(my.surv_hr_positive_exome_3 ~ surv_hr_positive_exome_3[,1])

surv_tnbc_exome_3 <- as.matrix(tnbc[,c("exome_cdr3_3", "d", "vital")])
my.surv_tnbc_exome_3 <- Surv(as.numeric(surv_tnbc_exome_3[,2]), as.numeric(surv_tnbc_exome_3[,3]))
my.fit_tnbc_exome_3 <- survfit(my.surv_tnbc_exome_3 ~ surv_tnbc_exome_3[,1])

plot(my.fit_her2_positive_exome_3, col=c(1,2,3), main="Survival of her2+ by exome cdr3", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low", "zero"), fill=c(1,2,3))
plot(my.fit_hr_positive_exome_3, col=c(1,2,3), main="Survival of hr+ by exome cdr3", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low", "zero"), fill=c(1,2,3))
plot(my.fit_tnbc_exome_3, col=c(1,2,3), main="Survival of tnbc by exome cdr3", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low", "zero"), fill=c(1,2,3))

# D&R by subtype

surv_her2_positive_exome_rna <- as.matrix(her2_positive[,c("exome_rna", "d", "vital")])
my.surv_her2_positive_exome_rna <- Surv(as.numeric(surv_her2_positive_exome_rna[,2]), as.numeric(surv_her2_positive_exome_rna[,3]))
my.fit_her2_positive_exome_rna <- survfit(my.surv_her2_positive_exome_rna ~ surv_her2_positive_exome_rna[,1])

surv_hr_positive_exome_rna <- as.matrix(hr_positive[,c("exome_rna", "d", "vital")])
my.surv_hr_positive_exome_rna <- Surv(as.numeric(surv_hr_positive_exome_rna[,2]), as.numeric(surv_hr_positive_exome_rna[,3]))
my.fit_hr_positive_exome_rna <- survfit(my.surv_hr_positive_exome_rna ~ surv_hr_positive_exome_rna[,1])

surv_tnbc_exome_rna <- as.matrix(tnbc[,c("exome_rna", "d", "vital")])
my.surv_tnbc_exome_rna <- Surv(as.numeric(surv_tnbc_exome_rna[,2]), as.numeric(surv_tnbc_exome_rna[,3]))
my.fit_tnbc_exome_rna <- survfit(my.surv_tnbc_exome_rna ~ surv_tnbc_exome_rna[,1])

plot(my.fit_her2_positive_exome_rna, col=c(1,2,3,4), main="Survival of her2+ by exome/rna presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("D-R-", "D-R+", "D+R-", "D+R+"), fill=c(1,2,3,4))
plot(my.fit_hr_positive_exome_rna, col=c(1,2,3,4), main="Survival of hr+ by exome/rna presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("D-R-", "D-R+", "D+R-", "D+R+"), fill=c(1,2,3,4))
plot(my.fit_tnbc_exome_rna, col=c(1,2,4), main="Survival of tnbc by exome/rna presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("D-R-", "D-R+", "D+R-", "D+R+"), fill=c(1,2,3,4))

# treg/cdr3 by subtype

surv_her2_positive_exome_treg <- as.matrix(her2_positive[,c("exome_treg", "d", "vital")])
my.surv_her2_positive_exome_treg <- Surv(as.numeric(surv_her2_positive_exome_treg[,2]), as.numeric(surv_her2_positive_exome_treg[,3]))
my.fit_her2_positive_exome_treg <- survfit(my.surv_her2_positive_exome_treg ~ surv_her2_positive_exome_treg[,1])

surv_hr_positive_exome_treg <- as.matrix(hr_positive[,c("exome_treg", "d", "vital")])
my.surv_hr_positive_exome_treg <- Surv(as.numeric(surv_hr_positive_exome_treg[,2]), as.numeric(surv_hr_positive_exome_treg[,3]))
my.fit_hr_positive_exome_treg <- survfit(my.surv_hr_positive_exome_treg ~ surv_hr_positive_exome_treg[,1])

surv_tnbc_exome_treg <- as.matrix(tnbc[,c("exome_treg", "d", "vital")])
my.surv_tnbc_exome_treg <- Surv(as.numeric(surv_tnbc_exome_treg[,2]), as.numeric(surv_tnbc_exome_treg[,3]))
my.fit_tnbc_exome_treg <- survfit(my.surv_tnbc_exome_treg ~ surv_tnbc_exome_treg[,1])

plot(my.fit_her2_positive_exome_treg, col=c(1,2,3,4), main="Survival of her2+ by exome/T-reg high/low", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("D+T+", "D-T+", "D+T-", "D-T-"), fill=c(1,2,3,4))
plot(my.fit_hr_positive_exome_treg, col=c(1,2,3,4), main="Survival of hr+ by exome/T-reg high/low", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("D+T+", "D-T+", "D+T-", "D-T-"), fill=c(1,2,3,4))
plot(my.fit_tnbc_exome_treg, col=c(1,2,3,4), main="Survival of tnbc by exome/T-reg high/low", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("D+T+", "D-T+", "D+T-", "D-T-"), fill=c(1,2,3,4))


survdiff(my.surv_her2_positive_exome_rna ~ surv_her2_positive_exome_rna[,1])
survdiff(my.surv_her2_positive_exome_treg ~ surv_her2_positive_exome_treg[,1])


# k-means clusters (only if ran cdr3_firehose_clustering before)
surv_clusters <- mydata[,c("fit.cluster", "days", "vital")]
surv_clusters$days <- surv_clusters$days/365

my.surv_clusters <- Surv(as.numeric(surv_clusters[,2]), as.numeric(surv_clusters[,3]))
my.fit_clusters <- survfit(my.surv_clusters ~ surv_clusters[,1])
plot(my.fit_clusters, col=c(1,2,3,4), main="Survival of non-her2+ by K-means clusters", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("1", "2", "3", "4"), fill=c(1,2,3,4))

survdiff(my.surv_clusters ~ surv_clusters[,1])

# hierarchical clusters
surv_clusters <- results_hc[,c("mycl", "days", "vital")]
surv_clusters$days <- surv_clusters$days/365

my.surv_clusters <- Surv(as.numeric(surv_clusters[,2]), as.numeric(surv_clusters[,3]))
my.fit_clusters <- survfit(my.surv_clusters ~ surv_clusters[,1])
plot(my.fit_clusters, col=c(1,2,3), main="Survival of non-her2+ by clusters", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("Cluster A", "Cluster B", "Cluster C"), fill=c(1,2,3))

survdiff(my.surv_clusters ~ surv_clusters[,1])

# poster plots
survdiff(my.surv_her2_positive_exome ~ surv_her2_positive_exome[,1])
#her2+
plot(my.fit_her2_positive_exome, col=c(1,3), main="HER2-positive by exome cdr3", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.6, c("high", "low"), fill=c(1,3))
legend(1,0.2, "p=0.0159", bty="n")