#set the table working on
library(dplyr)
library(ggplot2)

cdr3_results_current<-read.csv("cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_clustering_age_dedup_1078.txt", sep="\t")

#renaming groups
cdr3_results_current<-mutate(cdr3_results_current,clinical_group=ifelse(clinical_group=="her2-","HR+/Her2-",as.character(clinical_group)))
cdr3_results_current<-mutate(cdr3_results_current,clinical_group=ifelse(clinical_group=="her2+","Her2+",as.character(clinical_group)))
cdr3_results_current<-mutate(cdr3_results_current,clinical_group=ifelse(clinical_group=="triple-","TNBC",as.character(clinical_group)))

# calculating the iDNA score, rna_rpm, R and D groups
cdr3_results_current<-mutate(cdr3_results_current,exome_rpm=exome_imseq*10^6/exome_reads)
tmp<-select(filter(cdr3_results_current,exome_rpm!=0),patient_id,exome_rpm)
tmp<-mutate(tmp,iDNA_score=ntile(exome_rpm,10))
cdr3_results_current<-left_join(cdr3_results_current,tmp)

cdr3_results_current<-mutate(cdr3_results_current,iDNA_score=ifelse(is.na(iDNA_score),0,iDNA_score))
cdr3_results_current<-mutate(cdr3_results_current,rna_rpm=rna_imseq*10^6/rna_reads)
cdr3_results_current<-mutate(cdr3_results_current,D=ifelse(exome_rpm==0,"D-","D+"),R=ifelse(rna_rpm==0,"R-","R+"))
cdr3_results_current<-mutate(cdr3_results_current,blood_rpm=blood_cdr3*10^6/blood_reads)
cdr3_results_current<-mutate(cdr3_results_current,tcrb_rpm=tcrb_reads*10^6/rna_reads)

# cytolytic activity deciles
# calculating the iDNA score, rna_rpm, R and D groups
tmp<-select(filter(cdr3_results_current,cytolytic!=0),patient_id,cytolytic)
tmp<-mutate(tmp,cytolytic_score=ntile(cytolytic,10))
cdr3_results_current<-left_join(cdr3_results_current,tmp)
cdr3_results_current<-mutate(cdr3_results_current,cytolytic_score=ifelse(is.na(cytolytic_score),0,cytolytic_score))

#fraction of 5% TILs in each bin

cdr3_results_current %>% group_by(iDNA_score) %>% summarize(length(lymphocyte_percent[lymphocyte_percent>=5])/length(lymphocyte_percent))

#association with TIL content
wilcox.test(as.vector(unlist(select(filter(cdr3_results_current,lymphocyte_percent<=2),exome_rpm))),as.vector(unlist(select(filter(cdr3_results_current,lymphocyte_percent>2),exome_rpm))))

# comparing D+R+
#tmp<-select(mutate(cdr3_results_current,D=ifelse(exome_rpm==0,"D-","D+"),R=ifelse(rna_rpm==0,"R-","R+")),D,R)

table(interaction(cdr3_results_current$D,cdr3_results_current$R))

#################################################
#Plotting 
#######################################

#figure RNA vs DNA
ggplot(filter(cdr3_results_current,rna_rpm>0),aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),log10(rna_rpm)))+geom_boxplot(fill="grey")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("log10(RNA CDR3(RPM))")+theme_classic(base_size = 24)

# figure FracRNA pos vs DNA
frac_withRNA<-aggregate(rna_rpm~iDNA_score,data=filter(cdr3_results_current,!is.na(rna_rpm)),FUN=function(x) length(x[x>0.0255])/length(x))
ggplot(frac_withRNA,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors \n with high RNA CDR3 RPM")+theme_classic(base_size = 24)

#figure FracTIL vs DNA
frac_withTILs<-aggregate(lymphocyte_percent~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=5])/length(x))
ggplot(frac_withTILs,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),lymphocyte_percent))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥5% TILs")+theme_classic(base_size = 24)

# distribution by clinical groups
ggplot(filter(cdr3_results_current,clinical_group!="other"),aes(as.factor(iDNA_score),fill=clinical_group))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))

#figurre TIL ave vs DNA
#ggplot(cdr3_results_with_groups,aes(factor(reorder(tile.y,as.numeric(tile.y))),lymphocyte_percent))+geom_boxplot()+theme(text=element_text(size=32))+xlab("DNA CDR3(RPM) bin (increasing)")+ylab("percent TILs")

# figure FracPurity vs DNA
frac_pure<-aggregate(patient_purity~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=0.80])/length(x))
ggplot(frac_pure,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),patient_purity))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with ≥80% purity")+theme_classic(base_size = 24)

# figure FracNon vs DNA
frac_non<-aggregate(rate_non~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=8.323135e-07])/length(x))
ggplot(frac_non,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rate_non))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high nonsyn")+theme_classic(base_size = 24)

# distribution by mRNAseq clusters
cdr3_results_mRNA_clusters <- cdr3_results_current[!(is.na(cdr3_results_current$CLUS_mRNAseq_cNMF)),]
mode(cdr3_results_mRNA_clusters$CLUS_mRNAseq_cNMF) <- "character"
ggplot(cdr3_results_mRNA_clusters,aes(as.factor(iDNA_score),fill=CLUS_mRNAseq_cNMF))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))

# distribution by methylation clusters
cdr3_results_methylation_clusters <- cdr3_results_current[!(is.na(cdr3_results_current$CLUS_Methlyation_cNMF)),]
mode(cdr3_results_methylation_clusters$CLUS_Methlyation_cNMF) <- "character"
ggplot(cdr3_results_methylation_clusters,aes(as.factor(iDNA_score),fill=CLUS_Methlyation_cNMF))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("gray16", "gray32", "gray48", "red", "gray64", "gray80"))

# figure cytolytic vs DNA
frac_cyto<-aggregate(cytolytic~iDNA_score,data=cdr3_results_current,FUN=function(x) length(x[x>=72.48])/length(x))
ggplot(frac_cyto,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),cytolytic))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high cytolytic activity")+theme_classic(base_size = 24)

# distribution by gsva clusters
cdr3_results_current$gsva_cluster <- as.factor(cdr3_results_current$gsva_cluster)
ggplot(cdr3_results_current,aes(as.factor(iDNA_score),fill=gsva_cluster))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))

### associations uncovered by linear models ###
restriction_cdr3_columns <- c("iDNA_score", "foxp3", "rate_non", "rna_rpm", "blood_rpm", "clinical_group")
cdr3_results_cdr3_restricted <- cdr3_results_current[,restriction_cdr3_columns]
cdr3_results_cdr3_restricted <- cdr3_results_cdr3_restricted[complete.cases(cdr3_results_cdr3_restricted),]

### cdr3
# frac foxp3 vs idna
frac_foxp3<-aggregate(foxp3~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=53.21])/length(x))
ggplot(frac_foxp3,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),foxp3))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high foxp3")+theme_classic(base_size = 24)

# figure FracNon vs DNA
frac_non<-aggregate(rate_non~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=8.323135e-07])/length(x))
ggplot(frac_non,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rate_non))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high nonsyn")+theme_classic(base_size = 24)

# frac rna_rpm vs idna
frac_rna<-aggregate(rna_rpm~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=0.028650])/length(x))
ggplot(frac_rna,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high rna rpm")+theme_classic(base_size = 24)

# frac blood_rpm vs idna
frac_blood<-aggregate(blood_rpm~iDNA_score,data=cdr3_results_cdr3_restricted,FUN=function(x) length(x[x>=0.014250])/length(x))
ggplot(frac_blood,aes(factor(reorder(iDNA_score,as.numeric(iDNA_score))),blood_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("iDNA score")+ylab("Fraction of tumors with high blood_rpm")+theme_classic(base_size = 24)

# distribution by clinical groups
ggplot(filter(cdr3_results_cdr3_restricted,clinical_group!="other"),aes(as.factor(iDNA_score),fill=clinical_group))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("iDNA score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))

### cytolytic
restriction_cytolytic_columns <- c("cytolytic_score", "tcrb_rpm", "cd4", "cd8a", "foxp3", "rna_rpm", "clinical_group")
cdr3_results_cytolytic_restricted <- cdr3_results_current[,restriction_cytolytic_columns]
cdr3_results_cytolytic_restricted <- cdr3_results_cytolytic_restricted[complete.cases(cdr3_results_cytolytic_restricted),]

# frac tcrb vs cytolytic
frac_tcrb<-aggregate(tcrb_rpm~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=30.560])/length(x))
ggplot(frac_tcrb,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),tcrb_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high tcrb")+theme_classic(base_size = 24)

# figure cd4 vs cytolytic
frac_cd4<-aggregate(cd4~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=1102])/length(x))
ggplot(frac_cd4,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),cd4))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high cd4")+theme_classic(base_size = 24)

# frac cd8a vs cytolytic
frac_cd8a<-aggregate(cd8a~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=159.6])/length(x))
ggplot(frac_cd8a,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),cd8a))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high cd8a")+theme_classic(base_size = 24)

# frac foxp3 vs cytolytic
frac_foxp3<-aggregate(foxp3~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=56.3600])/length(x))
ggplot(frac_foxp3,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),foxp3))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high foxp3")+theme_classic(base_size = 24)

# frac rna_rpm vs cytolytic
frac_rna<-aggregate(rna_rpm~cytolytic_score,data=cdr3_results_cytolytic_restricted,FUN=function(x) length(x[x>=0.027990])/length(x))
ggplot(frac_rna,aes(factor(reorder(cytolytic_score,as.numeric(cytolytic_score))),rna_rpm))+geom_bar(stat="identity")+theme(text=element_text(size=32))+xlab("cytolytic score")+ylab("Fraction of tumors with high rna rpm")+theme_classic(base_size = 24)

# distribution by clinical groups
ggplot(filter(cdr3_results_cytolytic_restricted,clinical_group!="other"),aes(as.factor(cytolytic_score),fill=clinical_group))+geom_bar(position="fill")+theme_classic(base_size = 24)+ylab("Fraction of tumors")+xlab("cytolytic score")+theme(legend.title=element_blank())+scale_fill_manual(values=c("blue", "red", "black"))



# output supp table
cols <- c("sample_id", "exome_reads", "exome_imseq", "exome_clonotypes", "iDNA_score")
output <- cdr3_results_current[,cols]
output$sample_id <- substr(output$sample_id, 1, 16)
write.table(output, "tcrseq_brca_output_table.txt", quote=FALSE, sep="\t")