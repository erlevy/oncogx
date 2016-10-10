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

# gsva table
#gsva_table <- sapply(gsva_cluster_split, table)
#gsva_percent <- prop.table(gsva_table, 2)
#barplot(gsva_percent) #, col=c("darkblue", "red"), legend=c("ifng low", "ifng high"), ylab="Percent high/low")
#gsva_table_df <- melt(gsva_table)
#ggplot(gsva_table_df, aes(x=Var2, y=value, fill=Var1)) + geom_bar(stat="identity", position="fill") + labs(title="Fraction of patients by GSVA immune signature groupings", x="Cohort", y="Fraction") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# gsva facet plots
cancer_heatmap <- cancer_current[!is.na(cancer_current$T.regs),]
cancer_heatmap$gsva_annotated <- as.character(cancer_heatmap$gsva_cluster)
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="1","High",as.character(gsva_annotated)))
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="2","Low",as.character(gsva_annotated)))
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="3","Mixed_Innate",as.character(gsva_annotated)))
cancer_heatmap<-mutate(cancer_heatmap,gsva_annotated=ifelse(gsva_annotated=="4","Mixed_Adaptive",as.character(gsva_annotated)))

# to exclude: ACC, COAD, DLBC, ESCA, GBM, KICH, KIRC, KIRP, LGG, 
#             LIHC, MESO, PCPG, PRAD, READ, SARC, SKCM, TGCT, THCA, UCEC, USC
# most interesting ones: BLCA, BRCA, CESC, HNSC, LUAD, LUSC, OV, PAAD, STAD
cancer_heatmap$gsva_annotated <- factor(cancer_heatmap$gsva_annotated, levels=c("Low", "Mixed_Innate", "Mixed_Adaptive", "High"))
#lapply(split(cancer_heatmap, cancer_heatmap$cohort:cancer_heatmap$gsva_annotated), nrow)
cancer_heatmap <- filter(cancer_heatmap, cohort == "BLCA" | cohort == "BRCA" | cohort == "CESC" |
                           cohort == "HNSC" | cohort == "LUAD" | cohort == "LUSC" |
                           cohort == "OV" | cohort == "PAAD" | cohort == "STAD")

### FIGURE: gsva vs. cdr3 rpm facet plots ###
ggplot(cancer_heatmap, aes(x=gsva_annotated, y=log(exome_rpm,10))) + geom_boxplot() + labs(x="", y="Normalized tumor DNA CDR3 fraction (log10)") + facet_wrap(~cohort, ncol=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=14))
ggplot(cancer_heatmap, aes(x=gsva_annotated, y=log(rna_rpm,10))) + geom_boxplot() + labs(x="", y="Normalized tumor RNA CDR3 fraction (log10)") + facet_wrap(~cohort, ncol=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text=element_text(size=14))
#ggplot(cancer_heatmap, aes(x=gsva_annotated, y=log(blood_rpm,10))) + geom_boxplot() + labs(x="", y="Blood exome CDR3 RPM (log10)") + facet_wrap(~cohort, ncol=3) # + theme(text=element_text(size=16))


# gsva
cancer_heatmap <- cancer_current[!is.na(cancer_current$T.regs),]

signatures <- select(cancer_heatmap, cohort, T.regs:Neut)
cohorts <- signatures[,1]
signatures <- signatures[,-1]

signature_cols <- c("T Regs", "T CD8", "T CD4 Naive", "T CD4 Mem. Resting", "T CD4 Mem. Activated",
                    "T Folicular Helper", "T Gamma Delta", "B Naive", "B Memory", "Plasma", "NK Resting",
                    "NK Activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
                    "Dendritic Resting", "Dendritic Activated", "Mast Resting", "Mast Activated",
                    "Eosinophils", "Neutrophils")

# colors still not working, need to copy work form "brca_final_paper.R")
#cohort_colors <- color.scale(as.numeric(cohorts), color.spec="bluerange")
#gsva_colors <- color.scale(as.numeric(cancer_heatmap$gsva_cluster), color.spec="bluerange")
#cohort_colors <- matrix(c(cohort_colors, cohort_colors), ncol=2, byrow=FALSE)

heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none", Colv=NA, labRow="",
             labCol=signature_cols, RowSideColors=cohort_colors)

#heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none", Colv=NA, labRow="",
#             labCol=signature_cols, cexCol=2, margins=c(22,6))


# check OV trend
ov <- filter(cancer_heatmap, cohort=="OV")
ov_gsva <- split(ov, ov$gsva_annotated)
t.test(ov_gsva$Low$exome_rpm, ov_gsva$High$exome_rpm)

