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
library(PredictABEL)

results <- read.csv("/Users/Eric/tcrseq/new/final/pancan_results_8-28-2016.txt", sep="\t")
clonotypes <- read.csv("/Users/Eric/tcrseq/new/final/pancan_clonotypes_8-10-2016.txt", sep="\t")

# add in LUMP score
### purity
purity <- read.csv("/Users/Eric/tcrseq/new/purity/pancancer_purity_butte.txt", sep="\t")

# process purity matrix
purity <- as.matrix(purity)
rownames(purity) <- purity[,1]
purity <- purity[,-1]

# match samples
patient_purity <- c()
data_list <- c("LUMP")
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
lump_score <- patient_purity
results <- cbind(results, lump_score)
results$lump_score <- as.numeric(as.character(results$lump_score))

cdr3_results_current <- results
cancer_current <- cdr3_results_current

#cancer_current <- filter(cdr3_results_current, cohort != "DLBC")

#cancer_current <- filter(cdr3_results_current, cohort != "ACC", cohort != "CHOL",
#                         cohort != "DLBC", cohort != "KICH", cohort != "TGCT", cohort != "THYM")

cancer_current <- filter(cdr3_results_current, cohort != "THYM", cohort != "UVM", cohort != "CHOL",
                         cohort != "ACC", cohort != "CESC", cohort != "DLBC", cohort != "ESCA", 
                         cohort != "KICH", cohort != "LIHC", cohort != "MESO",
                         cohort != "PCPG", cohort != "TGCT", cohort != "UCS")

cancer_current$cohort <- factor(cancer_current$cohort)
cancer_current<-mutate(cancer_current,gsva_lowhigh=ifelse(gsva_cluster==1 | gsva_cluster==4,"1","0"))
cancer_current<-mutate(cancer_current,blood_level=ifelse(blood_rpm<median(cancer_current$blood_rpm, na.rm=TRUE),"0","1"))

tmp<-select(filter(cancer_current,blood_rpm!=0),patient_uuid,blood_rpm)
tmp<-mutate(tmp,iDNA_blood_score=ntile(blood_rpm,10))
cancer_current<-left_join(cancer_current,tmp)
cancer_current<-mutate(cancer_current,iDNA_blood_score=ifelse(is.na(iDNA_blood_score),0,iDNA_blood_score))

tmp<-select(filter(cancer_current,rate_non!=0),patient_uuid,rate_non)
tmp<-mutate(tmp,nonsyn_score=ntile(rate_non,10))
cancer_current<-left_join(cancer_current,tmp)
cancer_current<-mutate(cancer_current,nonsyn_score=ifelse(is.na(nonsyn_score),0,nonsyn_score))

# CDR3 nonzero logistic regression
cancer_current$gsva_cluster <- factor(cancer_current$gsva_cluster)
cancer_current <- within(cancer_current, cohort <- relevel(cohort, ref="LGG"))
mylogit <- glm(exome_level ~ cohort + days_to_birth + gender + gsva_cluster + patient_purity + rna_rpm + 
                 blood_rpm + tcrb_rpm + exome_reads, data = cancer_current, family = "binomial")

# CDR3 nonzero linear regression
mylm <- lm(exome_rpm ~ cohort + days_to_birth + gender + gsva_cluster + patient_purity + rna_rpm + 
                 blood_rpm + tcrb_rpm, data = cancer_current)

# GSVA lowhigh logistic regression
cancer_current$gsva_lowhigh <- factor(cancer_current$gsva_lowhigh)
cancer_current <- within(cancer_current, cohort <- relevel(cohort, ref="LGG"))
mylogit <- glm(gsva_lowhigh ~ cohort + days_to_birth + gender + blood_rpm,
               data = cancer_current, family = "binomial")

# Lump score linear regression
cancer_current$gsva_lowhigh <- factor(cancer_current$gsva_lowhigh)
cancer_current <- within(cancer_current, cohort <- relevel(cohort, ref="LGG"))
mylm <- lm(lump_score ~ cohort + days_to_birth + gender +blood_rpm,
           data = cancer_current)

### gsva low/high stepwise building
# start with non-cohort factors
cancer_current$gsva_lowhigh <- factor(cancer_current$gsva_lowhigh)
cancer_current <- within(cancer_current, cohort <- relevel(cohort, ref="LGG"))
cancer_current$age <- cancer_current$days_to_birth/-365
mylogit <- glm(gsva_lowhigh ~ age + gender + iDNA_blood_score,
               data = cancer_current, family = "binomial")
# age and gender significant
mylogit <- glm(gsva_lowhigh ~ age + gender + iDNA_blood_score + nonsyn_score,
               data = cancer_current, family = "binomial")
# age significant (with burden)
mylogit <- glm(gsva_lowhigh ~ age + gender + iDNA_blood_score + nonsyn_score + cohort,
               data = cancer_current, family = "binomial")
#write.table(round(exp(mylogit$coefficients),4), "block_3.txt", quote=FALSE, sep="\t")
# significant in age and burden, still none in blood
# pick top cohorts: BRCA, CESC, HNSC, KIRC, KIRP, LUAD, LUSC, OV, PRAD, STAD, TGCT
#cancer_reduced <- filter(cancer_current, cohort %in% c("BRCA", "CESC", "HNSC", "KIRC", "KIRP","LUAD", "LUSC", "OV", "PRAD", "STAD", "TGCT"))
cancer_reduced <- cancer_current
mylogit <- glm(gsva_lowhigh ~ age + gender + nonsyn_score + iDNA_blood_score*cohort,
               data = cancer_reduced, family = "binomial")
#write.table(round(exp(mylogit$coefficients),4), "block_4.txt", quote=FALSE, sep="\t")
# gender becomes significant, burden not, blood:CESC interaction significant (BRCA if only do interactions)
# also blood with LUSC and STAD

### lump score stepwise building
# start with non-cohort factors
cancer_current$gsva_lowhigh <- factor(cancer_current$gsva_lowhigh)
cancer_current <- within(cancer_current, cohort <- relevel(cohort, ref="LGG"))
mylm <- lm(lump_score ~ days_to_birth + gender + blood_rpm + rate_non,
               data = cancer_current)
# age, blood, and burden significant, not gender
mylm <- lm(lump_score ~ days_to_birth + gender + blood_rpm + rate_non + cohort,
           data = cancer_current)
# loses significance in age and blood, still in burden
# all cohorts significant except: KICH, PRAD
cancer_reduced <- filter(cancer_current, !(cohort %in% c("KICH", "PRAD")))
mylm <- lm(lump_score ~ days_to_birth + gender + rate_non + blood_rpm*cohort,
           data = cancer_reduced)
# burden still significant, all cohorts, LIHC:blood interaction


# gender becomes significant, blood:CESC interaction significant (BRCA if only do interactions)