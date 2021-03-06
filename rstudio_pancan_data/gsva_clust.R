library(dplyr)

results_all <- read.csv("/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-29-2016.txt", sep="\t")

results <- results_all[!is.na(results_all$T.regs),]
signatures <- select(results, T.regs:Neut)

# hierarchical cluster analysis
d <- dist(as.matrix(signatures))
hc <- hclust(d)
plot(hc, labels=FALSE)
mycl <- cutree(hc, h=max(hc$height/1.5))

results_hc <- cbind(results, mycl)

mycl_df <- mutate(results_hc,gsva_colors=ifelse(mycl==4,"red","black"))
cohort_colors <- cbind(mycl_df$gsva_colors, mycl_df$gsva_colors)
heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none", Colv=NA, labRow="", RowSideColors=cohort_colors)

# add clusters to the 1078 results table (6 without the expression data)
mycl_all <- c()
names(mycl) <- results_hc[,1]
for (i in 1:nrow(results_all))
{
  patient_uuid <- as.character(results_all$patient_uuid[i])
  mycl_all <- c(mycl_all, mycl[patient_uuid])
}
mycl_all <- as.matrix(as.numeric((mycl_all)))
rownames(mycl_all) <- results_all$patient_uuid
results_hc_all <- cbind(results_all, mycl_all)
colnames(results_hc_all)[44] <- "gsva_cluster"

write.table(results_hc_all, "/Users/Eric/tcrseq/new/processed/pancan_exome_rna_blood_tcrb_gsva_clinical_7-30-2016.txt", sep="\t", 
            row.names=FALSE, quote=FALSE)
