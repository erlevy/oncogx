library(dplyr)
library(ggplot2)

#copying the data into working table
lymphocyte_results_table_with_clonotypes <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_gsva_clinical_dedup.txt", sep="\t")

#pancancer<-lymphocyte_results_table
pancancer<-lymphocyte_results_table_with_clonotypes

# calculating RPM
pancancer<-mutate(pancancer,exome_rpm=exome_imseq*10^6/exome_reads)
ggplot(pancancer,aes(cancer_type,log10(exome_rpm),fill=cancer_type))+geom_violin()+geom_boxplot(width=0.1,fill="black")

#calculating iDNAscore

tmp<-select(filter(pancancer,exome_rpm!=0),patient_id,exome_rpm)
tmp<-mutate(tmp,iDNA_score=ntile(exome_rpm,10))
pancancer<-left_join(pancancer,tmp)
pancancer<-mutate(pancancer,iDNA_score=ifelse(is.na(iDNA_score),0,iDNA_score))

#aggregating
pancan_aggregate<-aggregate(patient_id~iDNA_score+cancer_type,data=pancancer,FUN=length)
cancer_type_agg<-aggregate(patient_id~cancer_type,data=pancancer,FUN=length)
cancer_type_agg<-rename(cancer_type_agg,Total=patient_id)
pancan_aggregate<-rename(pancan_aggregate,Num=patient_id)
pancan_aggregate<-left_join(pancan_aggregate,cancer_type_agg)
pancan_aggregate<-mutate(pancan_aggregate,fraction=Num/Total)

#aggregating CDR3>0
pancan_aggregate_wCDR3<-aggregate(patient_id~iDNA_score+cancer_type,data=filter(pancancer,iDNA_score>0),FUN=length)
cancer_type_agg_wCDR3<-aggregate(patient_id~cancer_type,data=filter(pancancer,iDNA_score>0),FUN=length)
cancer_type_agg_wCDR3<-rename(cancer_type_agg_wCDR3,Total=patient_id)
pancan_aggregate_wCDR3<-rename(pancan_aggregate_wCDR3,Num=patient_id)
pancan_aggregate_wCDR3<-left_join(pancan_aggregate_wCDR3,cancer_type_agg)
pancan_aggregate_wCDR3<-mutate(pancan_aggregate_wCDR3,fraction=Num/Total)

#plot CDR3 nonNull
tmp<-filter(pancan_aggregate,iDNA_score==0)
tmp$cancer_type<-factor(tmp$cancer_type, levels = tmp$cancer_type[order(tmp$fraction)])
tmp<-mutate(tmp,color=ifelse(cancer_type=="BRCA"|cancer_type=="PAAD"|cancer_type=="LUAD"|cancer_type=="PRAD","red","black"))
#tmp$color<-factor(tmp$color, levels = tmp$color[order(tmp$fraction)])
ggplot(tmp,aes(cancer_type,1-(Num/Total)))+geom_bar(stat="identity",fill=tmp$color[order(tmp$fraction)])+theme_classic(base_size=20)+theme(axis.line.x=element_line(size=2),axis.line.y=element_line(size=2),axis.ticks.x=element_line(size=2),axis.ticks.y=element_line(size=2),axis.text.x=element_text(angle=90, hjust = 1))+ylab("Fraction of patients with iDNA")

#plot the CDR3 positive
ggplot(filter(pancan_aggregate,iDNA_score==0),aes(cancer_type,Num/Total))+geom_bar(stat="identity")
ggplot(filter(pancan_aggregate,iDNA_score>0),aes(as.factor(iDNA_score),Num/Total))+geom_bar(stat="identity")+facet_wrap(~cancer_type)

#plot all
#ggplot(pancancer,aes(factor(iDNA_score),log10(exome_rpm),color=cancer_type))+geom_point()+geom_jitter()+theme_classic(base_size=20)+theme(axis.line.x=element_line(size=2),axis.line.y=element_line(size=2),axis.ticks.x=element_line(size=2),axis.ticks.y=element_line(size=2))+ylab("CDR3 abundance (log10(RPM))")+xlab("iDNA score")

#plot
tmp2<-aggregate(exome_reads~cancer_type,data=pancancer,FUN=mean)
cancer_type_agg$exome_reads<-tmp2$exome_reads
tmp3<-filter(pancan_aggregate,iDNA_score==0)
tmp3<-rename(tmp3,noCDR3=Num)
cancer_type_agg_iDNA0<-left_join(cancer_type_agg,tmp3)
ggplot(cancer_type_agg_iDNA0,aes(log10(exome_reads),1-fraction,label=cancer_type))+geom_point(aes(size=Total,color=tmp$color))+geom_text(hjust=-0.1,size=6)+theme_classic(base_size=32)+theme(axis.line.x=element_line(size=2),axis.line.y=element_line(size=2),axis.ticks.x=element_line(size=2),axis.ticks.y=element_line(size=2))+ylab("Fraction of patients with iDNA")+xlab("Number of reads (log10)")


#plot four
ggplot(filter(pancancer,(cancer_type=="PAAD" | cancer_type=="BRCA" | cancer_type=="LUAD" | cancer_type=="PRAD") & iDNA_score>0),aes(cancer_type,log10(exome_rpm),fill=cancer_type))+geom_violin()+geom_boxplot(width=0.1,fill="black")+theme_classic(base_size=20)+theme(axis.line.x=element_line(size=2),axis.line.y=element_line(size=2),axis.ticks.x=element_line(size=2),axis.ticks.y=element_line(size=2))+ylab("CDR3 abundance (log10(RPM))")+xlab("")

#plot diversity

ggplot(pancancer,aes(exome_clonotypes/exome_rpm,fill=cancer_type))+geom_boxplot()


