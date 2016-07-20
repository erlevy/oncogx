library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)
library(car)
library(gvlma)
library(MASS)
library(leaps)
library(relaimpo)
library(dplyr)
library(RColorBrewer)
library(plotrix)

gsva_dir <- "/Users/Eric/tcrseq/new/gsva"

files_gsva <- list.files(path=gsva_dir, pattern="*.txt", full.names=T, recursive=FALSE)

gsva1 <- read.csv(files_gsva[1], sep="\t")

heatmap.plus(t(as.matrix(gsva1)), col=bluered(51), scale="none", Colv=NA, labRow="")

heatmaps <- function(file_path)
{
  gsva_file <- read.csv(file_path, sep="\t")
  cancer_type <- basename(strsplit(file_path, '[.]')[[1]][1])
  samples <- paste("_n=", ncol(gsva_file), sep="")
  heatmap_title <- paste(cancer_type, samples, sep="")
#  heatmap.plus(t(as.matrix(gsva_file)), col=bluered(51), scale="none", Colv=NA, labRow="", main=heatmap_title)
  hist(as.matrix(gsva_file[1,]))
}

lapply(files_gsva, heatmaps)


results <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_gsva_clinical_dedup.txt", sep="\t")
results_heatmap <- results[!is.na(results$T.regs),]

signatures <- results_heatmap[,c(3,12:33)]
signatures <- unique(signatures)
cohorts <- signatures[,1]
signatures <- signatures[,-1]

signature_cols <- c("T Regs", "T CD8", "T CD4 Naive", "T CD4 Mem. Resting", "T CD4 Mem. Activated",
                    "T Folicular Helper", "T Gamma Delta", "B Naive", "B Memory", "Plasma", "NK Resting",
                    "NK Activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
                    "Dendritic Resting", "Dendritic Activated", "Mast Resting", "Mast Activated",
                    "Eosinophils", "Neutrophils")

cohort_numbers <- as.numeric(cohorts)
#cohort_colors <- color.scale(cohort_numbers, c(0,1,1),c(1,1,0),0)

cohort_colors <- cbind(cohort_colors, cohort_colors)

heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none", Colv=NA, labRow="", labCol=signature_cols, cexCol=2, margins=c(22,6))
             
heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=cohort_colors,
             scale="none", Colv=NA, labRow="", labCol=signature_cols)

heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none",
             Colv=NA, labRow="", cexCol=4.5, margins=c(36,6), labCol=signature_cols)



# # heatmap of gsva clusters
# #results_clustering <- results[!is.na(results$gsva_cluster),]
# #results_clustering <- results_clustering[!is.na(results_clustering$rna_reads),]
# #results_clustering <- results_clustering[which(results_clustering$clinical_group!="other"),]
# #results_clustering <- results_clustering[!is.na(results_clustering$clinical_group),]
# #col_mat <- matrix(nrow=nrow(results_clustering), ncol=4)
# #rna_freq <- results_clustering$rna_imseq/results_clustering$rna_reads
# #rna_med <- median(rna_freq)
# for (i in 1:nrow(col_mat))
# {
#   group <- as.character(results_clustering$clinical_group[i])
#   exome <- results_clustering$exome_imseq[i]
#   rna <- rna_freq[i]
#   cluster_patient <- results_clustering$gsva_cluster[i]
#   if (exome==0) {col_mat[i,2]="white"}
#   else {col_mat[i,2]="black"}
#   if (rna<=rna_med) {col_mat[i,3]="white"}
#   else {col_mat[i,3]="purple"}
#   if (group=="her2+") {col_mat[i,1]="red"}
#   else if (group=="her2-") {col_mat[i,1]="green"}
#   else if (group=="triple-") {col_mat[i,1]="blue"}
#   if (cluster_patient==1) {col_mat[i,4]="red"}
#   else if (cluster_patient==2) {col_mat[i,4]="green"}
#   else if (cluster_patient==3) {col_mat[i,4]="blue"}
# }
# 
# ggplotColours <- function(n=6, h=c(0, 360) +15){
#   if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
#   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }
# 
# brewer.pal(4, "YlOrRd")
# 
# for (i in 1:nrow(col_mat))
# {
#   exome <- results_clustering$exome_imseq[i]
#   cluster_patient <- results_clustering$gsva_cluster[i]
#   #  cluster_patient <- results_clustering$mycl[i]
#   if (exome==0) {col_mat[i,2]="white"}
#   else {col_mat[i,2]="black"}
#   if (cluster_patient==1) {col_mat[i,4]="#FECC5C"}  # middle/innate
#   else if (cluster_patient==2) {col_mat[i,4]="#FD8D3C"} # middle/adaptive
#   else if (cluster_patient==3) {col_mat[i,4]="#E31A1C"} # high
#   else if (cluster_patient==4) {col_mat[i,4]="#FFFFB2"}  # low
# }
# 
# col_mat[,2] <- col_mat[,4]
# col_mat[,3] <- col_mat[,4]
# col_mat[1072,1] <- NA
# 
# signatures <- results_clustering[,50:71] # new cdr3_results with dedup
# rownames(signatures) <- results_clustering[,1]
# signatures <- scale(signatures)
# 
# signature_cols <- c("T Regs", "T CD8", "T CD4 Naive", "T CD4 Mem. Resting", "T CD4 Mem. Activated",
#                     "T Folicular Helper", "T Gamma Delta", "B Naive", "B Memory", "Plasma", "NK Resting",
#                     "NK Activated", "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
#                     "Dendritic Resting", "Dendritic Activated", "Mast Resting", "Mast Activated",
#                     "Eosinophils", "Neutrophils")
# 
# heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none",
#              Colv=NA, labRow="", cexCol=4.5, margins=c(36,6), labCol=signature_cols)
# 

# es.dif_col <- c()
# es.dif_all <- gsub(".", "-", colnames(es.dif), fixed=TRUE)
# es.dif_all <- substr(es.dif_all,1,16)
# #es.dif_all <- gsub("X", "", es.dif_all, fixed=TRUE)
# for (i in 1:nrow(cdr3))
# {
#   sample_id <- substr(as.character(cdr3$sample_id[i]),1,16)
#   es.dif_index <- which(sample_id==es.dif_all)
#   if (length(es.dif_index)>0)
#   {
#     es.dif_value <- es.dif[,es.dif_index]
#     es.dif_col <- rbind(es.dif_col, es.dif_value)
#   }
#   else
#   {
#     es.dif_col <- rbind(es.dif_col, rep(NA, 22))
#   }
# }
# rownames(es.dif_col) <- seq(1,nrow(cdr3))
# cdr3_new <- cbind(cdr3, es.dif_col)
# write.csv(cdr3_new, "cdr3_results_with_expression_blood_ptprc_burden_purity_gsva_1078.txt", quote=FALSE, row.names=FALSE)
