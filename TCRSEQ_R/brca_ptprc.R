
#results <- read.csv("cdr3_results_with_groups_complete.txt", sep="\t")
results <- read.csv("cdr3_results_with_expression_blood_1078.txt", sep="\t")
expression <- read.csv("BRCA_rnaseqv2_RSEM_isoforms_normalized_PTPRC.txt", sep="\t")
ucsc <- read.csv("ucsc_genes_ptprc.txt", sep="\t")

expression <- as.matrix(expression)
rownames(expression) <- expression[,1]
expression <- expression[,-1]
mode(expression) <- 'numeric'
expression <- t(expression)
rownames(expression) <- gsub("X", "", fixed=TRUE, rownames(expression))
rownames(expression) <- gsub(".", "-", fixed=TRUE, rownames(expression))
rownames(expression) <- substr(rownames(expression), start=1, stop=15)

patient_expression <- c()
for (i in 1:nrow(results))
{
  patient <- as.character(substr(results$sample_id[i],1,15))
  matching <- which(rownames(expression)==patient)
  if (length(matching)>0)
  {
    matching_expression <- expression[matching,]
  }
  else
  {
    matching_expression <- rep(NA,10)
  }
  patient_expression <- rbind(patient_expression, matching_expression)
}
results <- cbind(results, patient_expression)

results_heatmap <- results[,c("patient_id", "exome_imseq", "clinical_group", "uc001guq.2", "uc001gur.1",       
                              "uc001gus.1", "uc001gut.1", "uc001guu.1", "uc001guv.1", "uc001guw.1",         
                              "uc009wze.1", "uc009wzf.1", "uc010ppg.1")]
results_heatmap <- results_heatmap[complete.cases(results_heatmap),]
col_mat <- matrix(nrow=nrow(results_heatmap), ncol=2)
for (i in 1:nrow(col_mat))
{
  group <- as.character(results_heatmap$clinical_group[i])
  exome <- as.numeric(results_heatmap$exome_imseq[i])
  if (exome==0) {col_mat[i,1]="black"}
  else {col_mat[i,1]="red"}
  if (is.na(group)) {col_mat[i,2]="grey"}
  else if (group=="her2+") {col_mat[i,2]="red"}
  else if (group=="her2-") {col_mat[i,2]="green"}
  else if (group=="tnbc") {col_mat[i,2]="blue"}
  else {col_mat[i,2]="black"}
}

signatures <- results_heatmap[,-c(1,2,3)]
rownames(signatures) <- results_heatmap[,1]
heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="row")
heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none")
heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="col")

# rename isoform ids
ids <- c("Q0VAE8", "P08575", "uc001gus.1", " P08575-2", "F5GZM5", "uc001guv.1", "uc001guw.1", "Q5T9M4",
         "E9PKH0", "F5GXZ3")
colnames(results)[28:37] <- ids
write.table(results, "cdr3_results_with_expression_blood_ptprc_1078.txt", quote=FALSE, sep="\t", row.names=FALSE)
