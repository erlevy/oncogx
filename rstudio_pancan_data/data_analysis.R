library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)
library(reshape)

results <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-30-2016.txt", sep="\t")

### FROM OLD BRCA ###
#set the table working on
library(dplyr)
library(ggplot2)

cdr3_results_current<-results

#renaming groups
#cdr3_results_current<-mutate(cdr3_results_current,clinical_group=ifelse(clinical_group=="her2-","HR+/Her2-",as.character(clinical_group)))
#cdr3_results_current<-mutate(cdr3_results_current,clinical_group=ifelse(clinical_group=="her2+","Her2+",as.character(clinical_group)))
#cdr3_results_current<-mutate(cdr3_results_current,clinical_group=ifelse(clinical_group=="triple-","TNBC",as.character(clinical_group)))

# calculating the iDNA score, rna_rpm, R and D groups
cdr3_results_current<-mutate(cdr3_results_current,exome_rpm=exome_cdr3*10^6/exome_reads)
tmp<-select(filter(cdr3_results_current,exome_rpm!=0),patient_uuid,exome_rpm)
tmp<-mutate(tmp,iDNA_score=ntile(exome_rpm,10))
cdr3_results_current<-left_join(cdr3_results_current,tmp)
cdr3_results_current<-mutate(cdr3_results_current,iDNA_score=ifelse(is.na(iDNA_score),0,iDNA_score))

cdr3_results_current<-mutate(cdr3_results_current,rna_rpm=rna_cdr3*10^6/rna_reads)
cdr3_results_current<-mutate(cdr3_results_current,D=ifelse(exome_rpm==0,"D-","D+"),R=ifelse(rna_rpm==0,"R-","R+"))
#cdr3_results_current<-mutate(cdr3_results_current,blood_rpm=blood_cdr3*10^6/blood_reads)
#cdr3_results_current<-mutate(cdr3_results_current,tcrb_rpm=tcrb_reads*10^6/rna_reads)

#fraction of 5% TILs in each bin

#cdr3_results_current %>% group_by(iDNA_score) %>% summarize(length(percent_tils[percent_tils>=5])/length(percent_tils))

#association with TIL content
#wilcox.test(as.vector(unlist(select(filter(cdr3_results_current,percent_tils<=2),exome_rpm))),as.vector(unlist(select(filter(cdr3_results_current,percent_tils>2),exome_rpm))))

# comparing D+R+
#tmp<-select(mutate(cdr3_results_current,D=ifelse(exome_rpm==0,"D-","D+"),R=ifelse(rna_rpm==0,"R-","R+")),D,R)

#table(interaction(cdr3_results_current$D,cdr3_results_current$R))

#################################################
#Plotting 
#######################################

#figure RNA vs DNA
ggplot(filter(cdr3_results_current,rna_rpm>0),aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),log10(rna_rpm)))+geom_boxplot(fill="grey")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("log10(RNA CDR3(RPM))")+theme_classic(base_size = 24)

# figure FracRNA pos vs DNA
frac_withRNA<-aggregate(rna_rpm~iDNA_score,data=filter(cdr3_results_current,!is.na(rna_rpm)),FUN=function(x) length(x[x>0.0255])/length(x))
ggplot(frac_withRNA,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors \n with high RNA CDR3 RPM")+theme_classic(base_size = 24)

#figure FracTIL vs DNA
frac_withTILs<-aggregate(percent_tils~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=5])/length(x))
ggplot(frac_withTILs,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),percent_tils))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥5% TILs")+theme_classic(base_size = 24)

# distribution by clinical groups
ggplot(cdr3_results_current,aes(as.factor(iDNA_score),fill=cohort))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())

#figurre TIL ave vs DNA
#ggplot(cdr3_results_with_groups,aes(factor(reorder(tile.y,as.numeric(tile.y))),percent_tils))+geom_boxplot()+theme(text=element_text(size=32))+xlab("DNA CDR3(RPM) bin (increasing)")+ylab("percent TILs")

# # figure FracPurity vs DNA
# frac_pure<-aggregate(patient_purity~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=0.80])/length(x))
# ggplot(frac_pure,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),patient_purity))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥80% purity")+theme_classic(base_size = 24)
# 
# # figure FracNon vs DNA
# frac_non<-aggregate(rate_non~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=8.323135e-07])/length(x))
# ggplot(frac_non,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rate_non))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high nonsyn")+theme_classic(base_size = 24)
# 
# # distribution by mRNAseq clusters
# cdr3_results_mRNA_clusters <- cdr3_results_current[!(is.na(cdr3_results_current$CLUS_mRNAseq_cNMF)),]
# mode(cdr3_results_mRNA_clusters$CLUS_mRNAseq_cNMF) <- "character"
# ggplot(cdr3_results_mRNA_clusters,aes(as.factor(iDNA_score),fill=CLUS_mRNAseq_cNMF))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))
# 
# # distribution by methylation clusters
# cdr3_results_methylation_clusters <- cdr3_results_current[!(is.na(cdr3_results_current$CLUS_Methlyation_cNMF)),]
# mode(cdr3_results_methylation_clusters$CLUS_Methlyation_cNMF) <- "character"
# ggplot(cdr3_results_methylation_clusters,aes(as.factor(iDNA_score),fill=CLUS_Methlyation_cNMF))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("gray16", "gray32", "gray48", "red", "gray64", "gray80"))
# 
# # figure cytolytic vs DNA
# frac_cyto<-aggregate(cytolytic~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=72.48])/length(x))
# ggplot(frac_cyto,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),cytolytic))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high cytolytic activity")+theme_classic(base_size = 24)
# 
# # distribution by gsva clusters
# cdr3_results_current$gsva_cluster <- as.factor(cdr3_results_current$gsva_cluster)
# ggplot(cdr3_results_current,aes(as.factor(iDNA_score),fill=gsva_cluster))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))
# 
# ### associations uncovered by linear models ###
# restriction_cdr3_columns <- c("iDNA_score", "foxp3", "rate_non", "rna_rpm", "blood_rpm", "clinical_group")
# cdr3_results_cdr3_restricted <- cdr3_results_current[,restriction_cdr3_columns]
# cdr3_results_cdr3_restricted <- cdr3_results_cdr3_restricted[complete.cases(cdr3_results_cdr3_restricted),]
# 
# ### cdr3
# # frac foxp3 vs idna
# frac_foxp3<-aggregate(foxp3~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=53.21])/length(x))
# ggplot(frac_foxp3,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),foxp3))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high foxp3")+theme_classic(base_size = 24)
# 
# # figure FracNon vs DNA
# frac_non<-aggregate(rate_non~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=8.323135e-07])/length(x))
# ggplot(frac_non,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rate_non))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high nonsyn")+theme_classic(base_size = 24)
# 
# # frac rna_rpm vs idna
# frac_rna<-aggregate(rna_rpm~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=0.028650])/length(x))
# ggplot(frac_rna,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high rna rpm")+theme_classic(base_size = 24)
# 
# # frac blood_rpm vs idna
# frac_blood<-aggregate(blood_rpm~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=0.014250])/length(x))
# ggplot(frac_blood,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),blood_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high blood_rpm")+theme_classic(base_size = 24)
# 
# # distribution by clinical groups
# ggplot(filter(cdr3_results_cdr3_restricted,clinical_group!="other"),aes(as.factor(iDNA_score),fill=clinical_group))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))
# 
# ### cytolytic
# restriction_cytolytic_columns <- c("cytolytic_score", "tcrb_rpm", "cd4", "cd8a", "foxp3", "rna_rpm", "clinical_group")
# cdr3_results_cytolytic_restricted <- cdr3_results_current[,restriction_cytolytic_columns]
# cdr3_results_cytolytic_restricted <- cdr3_results_cytolytic_restricted[complete.cases(cdr3_results_cytolytic_restricted),]
# 
# # frac tcrb vs cytolytic
# frac_tcrb<-aggregate(tcrb_rpm~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=30.560])/length(x))
# ggplot(frac_tcrb,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),tcrb_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high tcrb")+theme_classic(base_size = 24)
# 
# # figure cd4 vs cytolytic
# frac_cd4<-aggregate(cd4~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=1102])/length(x))
# ggplot(frac_cd4,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),cd4))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high cd4")+theme_classic(base_size = 24)
# 
# # frac cd8a vs cytolytic
# frac_cd8a<-aggregate(cd8a~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=159.6])/length(x))
# ggplot(frac_cd8a,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),cd8a))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high cd8a")+theme_classic(base_size = 24)
# 
# # frac foxp3 vs cytolytic
# frac_foxp3<-aggregate(foxp3~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=56.3600])/length(x))
# ggplot(frac_foxp3,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),foxp3))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high foxp3")+theme_classic(base_size = 24)
# 
# # frac rna_rpm vs cytolytic
# frac_rna<-aggregate(rna_rpm~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=0.027990])/length(x))
# ggplot(frac_rna,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high rna rpm")+theme_classic(base_size = 24)
# 
# # distribution by clinical groups
# ggplot(filter(cdr3_results_cytolytic_restricted,clinical_group!="other"),aes(as.factor(cytolytic_score),fill=clinical_group))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("cytolytic score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))
# 
# 
# # output supp table
# cols <- c("sample_id", "exome_reads", "exome_imseq", "exome_clonotypes", "iDNA_score")
# output <- cdr3_results_current[,cols]
# output$sample_id <- substr(output$sample_id, 1, 16)
# write.table(output, "tcrseq_brca_output_table.txt", quote=FALSE, sep="\t")

### intial analysis
table(cdr3_results_current$cohort)
sum(!is.na(cdr3_results_current$exome_cdr3))
sum(!is.na(cdr3_results_current$rna_cdr3))
summary(cdr3_results_current$exome_cdr3)
summary(cdr3_results_current$rna_cdr3)
length(which(cdr3_results_current$exome_cdr3==0))
length(which(cdr3_results_current$rna_cdr3==0))
table(cdr3_results_current[,c("D", "R")])

cohorts <- split(cdr3_results_current, cdr3_results_current$cohort)
cohorts_limited_df <- rbind(cohorts$BLCA, cohorts$BRCA, cohorts$COAD, cohorts$GBM,
                        cohorts$HNSC, cohorts$LGG, cohorts$LUAD, cohorts$LUSC,
                        cohorts$OV, cohorts$PRAD, cohorts$STAD, cohorts$THCA)
cohorts_limited_list <- list(cohorts$BLCA, cohorts$BRCA, cohorts$COAD, cohorts$GBM,
                            cohorts$HNSC, cohorts$LGG, cohorts$LUAD, cohorts$LUSC,
                            cohorts$OV, cohorts$PRAD, cohorts$STAD, cohorts$THCA)
cohorts_limited_df$cohort <- factor(cohorts_limited_df$cohort)
boxplot(exome_reads ~ cohort, data=cdr3_results_current)
boxplot(exome_cdr3 ~ cohort, data=cdr3_results_current)
boxplot(log(rna_cdr3, 10) ~ cohort, data=cdr3_results_current)
boxplot(log(rna_rpm, 10) ~ cohort, data=cdr3_results_current)
boxplot(percent_tils ~ cohort, data=cdr3_results_current)

cohorts_limited_df$exome_cdr3[which(cohorts_limited$exome_cdr3==0)] <- NA
cohorts_limited_df$rna_cdr3[which(cohorts_limited$rna_cdr3==0)] <- NA
boxplot(log(exome_cdr3,10) ~ cohort, data=cohorts_limited_df)
boxplot(log(rna_cdr3,10) ~ cohort, data=cohorts_limited_df)

cdr3_profile <- function(cdr3_results)
{
  c(mean(cdr3_results$exome_cdr3, na.rm = TRUE), max(cdr3_results$exome_cdr3, na.rm = TRUE),
    mean(cdr3_results$rna_cdr3, na.rm = TRUE), max(cdr3_results$rna_cdr3, na.rm = TRUE),
    length(which(cdr3_results$exome_cdr3==0)), length(which(cdr3_results$rna_cdr3==0)))
}

cohorts_limited_profiling <- lapply(cohorts_limited_list, cdr3_profile)
df <- do.call("rbind", cohorts_limited_profiling)
write.table(df, "12_cohorts_profiling.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

### survival
### process days
cdr3_results_surv <- c()
for (i in 1:nrow(cdr3_results_current))
{
  surv_row <- rep(NA, 3)
  idna <- cdr3_results_current$iDNA_score[i]
  days_to_death <- cdr3_results_current$days_to_death[i]
  days_to_last_followup <- cdr3_results_current$days_to_last_followup[i]
  vital <- NA
  days <- NA
  if (idna>0) {idna <- 1}
  if (!is.na(days_to_death)) {
    days <- days_to_death
    vital <- 1
  } else if (!is.na(days_to_last_followup)) {
    days <- days_to_last_followup
    vital <- 0
    }
  surv_row <- c(idna, days, vital)
  cdr3_results_surv <- rbind(cdr3_results_surv, surv_row)
}
colnames(cdr3_results_surv) <- c("iDNA_nonzero", "days_compound", "vital_status")
rownames(cdr3_results_surv) <- seq(1:nrow(cdr3_results_surv))
cdr3_results_surv <- cbind(cdr3_results_current, cdr3_results_surv)
results_out_surv <- cbind(results, cdr3_results_surv$days_compound, cdr3_results_surv$vital_status)
colnames(results_out_surv)[45] <- "days_compound"
colnames(results_out_surv)[46] <- "vital_status"

write.table(results_out_surv, "/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_surv_7-30-2016.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)


surv_results <- function(surv_matrix)
{
  surv <- surv_matrix[,c("iDNA_nonzero", "days_compound", "vital_status")]
  surv$days_compound <- surv_matrix$days_compound/365
  my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
  my.fit <- survfit(my.surv ~ surv[,1])
#  plot(my.fit, col=c(1,2), main="Survival of HR+ by CDR3 groups", xlab="Survival Time (years)", ylab="Survival Probability")
#  legend(0.5,0.4, c("nonzero", "zero"), fill=c(1,2), cex=0.75)
  coxph(my.surv ~ surv[,1], method="breslow")  
}

types <- split(cdr3_results_surv, cdr3_results_surv$cohort)
types_surv <- lapply(types, surv_results)

# looks like only COAD is significant
# close: SKCM, CESC, ACC
surv_matrix <- types$COAD
surv <- surv_matrix[,c("iDNA_nonzero", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of COAD by CDR3 presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(5,0.4, c("zero", "nonzero"), fill=c(1,2), cex=0.75)
coxph(my.surv ~ surv[,1], method="breslow")  

surv_matrix <- types$SKCM
surv <- surv_matrix[,c("iDNA_nonzero", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of SKCM by CDR3 presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(1,0.4, c("zero", "nonzero"), fill=c(1,2), cex=0.75)
coxph(my.surv ~ surv[,1], method="breslow")  

surv_matrix <- types$CESC
surv <- surv_matrix[,c("iDNA_nonzero", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of CESC by CDR3 presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(8,0.4, c("zero", "nonzero"), fill=c(1,2), cex=0.75)
coxph(my.surv ~ surv[,1], method="breslow")  

surv_matrix <- types$ACC
surv <- surv_matrix[,c("iDNA_nonzero", "days_compound", "vital_status")]
surv$days_compound <- surv_matrix$days_compound/365
my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of ACC by CDR3 presence", xlab="Survival Time (years)", ylab="Survival Probability")
legend(5,0.4, c("zero", "nonzero"), fill=c(1,2), cex=0.75)
coxph(my.surv ~ surv[,1], method="breslow")  

### gsva pca
cdr3_results_gsva <- cdr3_results_current[,c(3, 16:37)]
cdr3_results_gsva <- cdr3_results_gsva[complete.cases(cdr3_results_gsva),]
cdr3_results_gsva_cohorts <- cdr3_results_gsva$cohort
cdr3_results_gsva <- t(cdr3_results_gsva[,-1])
pca <- prcomp(cdr3_results_gsva)
plot(pca[[2]][,1], pca[[2]][,2], col=cdr3_results_gsva_cohorts)

### rna vs. dna clonal diversity
cdr3_results_clonality <-mutate(cdr3_results_current,exome_clonality=exome_clonotypes/exome_rpm)
cdr3_results_clonality <-mutate(cdr3_results_clonality,rna_clonality=rna_clonotypes/rna_rpm)
plot(cdr3_results_clonality$exome_rpm, cdr3_results_clonality$exome_clonotypes)
plot(cdr3_results_clonality$rna_rpm, cdr3_results_clonality$rna_clonotypes)
plot(cdr3_results_clonality$exome_clonality, cdr3_results_clonality$rna_clonality)

cohorts_limited_df <- rbind(cohorts$BRCA, cohorts$HNSC,cohorts$LUAD, cohorts$LUSC,
                            cohorts$OV, cohorts$STAD)
cohorts_limited_df$cohort <- factor(cohorts_limited_df$cohort)
cdr3_results_clonality <-mutate(cohorts_limited_df,exome_clonality=exome_clonotypes/exome_rpm)
cdr3_results_clonality <-mutate(cdr3_results_clonality,rna_clonality=rna_clonotypes/rna_rpm)
plot(cdr3_results_clonality$exome_clonality, cdr3_results_clonality$rna_clonality,
     col=c("blue", "red", "black", "purple", "green", "lightblue"), pch=16,
     main="Clonal diversity of RNA vs. DNA", xlab="Exome clonal diversity", ylab="RNA clonal diversity")
legend(x=400, y=500, legend=levels(cdr3_results_clonality$cohort), pch=16,
       col=c("blue", "red", "black", "purple", "green", "lightblue"))


# add purity

### purity
purity <- read.csv("/Users/Eric/tcrseq/new/purity/pancancer_purity_butte.txt", sep="\t")

# process purity matrix
purity <- as.matrix(purity)
rownames(purity) <- purity[,1]
purity <- purity[,-1]

# match samples
patient_purity <- c()
data_list <- c("CPE")
for (i in 1:nrow(results))
{
  sample_id <- substr(as.character(results$sample_barcode[i]),1,16)
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

write.table(results_purity, "/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_surv_8-3-2016.txt", quote=FALSE, sep="\t", row.names=FALSE)


