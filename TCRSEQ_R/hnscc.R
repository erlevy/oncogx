library(OIsurv)
library(cgdsr)
library(GSVA)
library(heatmap.plus)
library(gplots)

#expression <- read.csv("/Users/Eric/tcga/hnscc/GSE40774_series_matrix_processed.csv", row.names=1)
key <- read.csv("/Users/Eric/hnscc/geo/agilent_ids_no_excel.txt", sep="\t", row.names=1)
#expression <- read.csv("BRCA_expression_firehose_patient_ids.csv")
#cdr3 <- read.csv("cdr3_results_with_groups.txt", sep="\t")

# process expression matrix
# expression_new <- c()
# expression_rows <- rownames(expression)
# for (i in 1:nrow(expression))
# {
#   expression_row <- expression[i,]
#   agilent_id <- expression_rows[i]
#   matching <- as.character(key[agilent_id,]$GENE_SYMBOL)
#   if (matching!="")
#   {
#     rownames(expression_row) <- matching
#     expression_new <- rbind(expression_new, expression_row)
#   }
# }
#expression_new <- read.csv("/Users/Eric/hnscc/geo/hnscc_expression_gene_ids.csv", row.names=1)
expression_new <- read.csv("/Users/Eric/hnscc/FromUC/UC-HNSC-Agilent/ProcessedData.gene.v2.134.txt", sep="\t", row.names=1)

# get cibersort list of gene signatures
cibersort <- read.csv("cibersort_signatures.txt", sep="\t")
tregs <- cibersort$genes[which(cibersort$T.cells.regulatory..Tregs.==1)]
tcd8 <- cibersort$genes[which(cibersort$T.cells.CD8==1)]
tcd4n <- cibersort$genes[which(cibersort$T.cells.CD4.naive==1)]
tcd4mr <- cibersort$genes[which(cibersort$T.cells.CD4.memory.resting==1)]
tcd4ma <- cibersort$genes[which(cibersort$T.cells.CD4.memory.activated==1)]
tcfh <- cibersort$genes[which(cibersort$T.cells.follicular.helper==1)]
tcgd <- cibersort$genes[which(cibersort$T.cells.gamma.delta==1)]

bcn <- cibersort$genes[which(cibersort$B.cells.naive==1)]
bcm <- cibersort$genes[which(cibersort$B.cells.memory==1)]
pc <- cibersort$genes[which(cibersort$Plasma.cells==1)]
nkr <- cibersort$genes[which(cibersort$NK.cells.resting==1)]
nka <- cibersort$genes[which(cibersort$NK.cells.activated==1)]
mono <- cibersort$genes[which(cibersort$Monocytes==1)]
macro0 <- cibersort$genes[which(cibersort$Macrophages.M0==1)]
macro1 <- cibersort$genes[which(cibersort$Macrophages.M1==1)]
macro2 <- cibersort$genes[which(cibersort$Macrophages.M2==1)]
denr <- cibersort$genes[which(cibersort$Dendritic.cells.resting==1)]
dena <- cibersort$genes[which(cibersort$Dendritic.cells.activated==1)]
mastr <- cibersort$genes[which(cibersort$Mast.cells.resting==1)]
masta <- cibersort$genes[which(cibersort$Mast.cells.activated==1)]
eos <- cibersort$genes[which(cibersort$Eosinophils==1)]
neut <- cibersort$genes[which(cibersort$Neutrophils==1)]
all_signatures <- list(tregs, tcd8, tcd4n, tcd4mr, tcd4ma, tcfh, tcgd, bcn, bcm, pc, nkr,
                       nka, mono, macro0, macro1, macro2, denr, dena, mastr, masta, eos, neut)
all_signature_names <- c("T regs", "T-CD8d", "T-CD4n", "T-CD4mr", "T-CD4ma", "T-hf", "T-gd", "B-n", "B-m", "Plasma",
                         "NK-r", "NK-a", "Mono", "Macro-0", "Macro-1", "Macro-2", "Den-r", "Den-a", "Mast-r", "Mast-a",
                         "Eos", "Neut")
names(all_signatures) <- all_signature_names

# process expression matrix
expression <- expression_new
expression <- as.matrix(expression)
mode(expression) <- 'numeric'
gene_names <- rownames(expression)
all_signatures <- lapply(all_signatures, as.character)

# run gsva
es.max_full <- gsva(expression, all_signatures, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif_full <- gsva(expression, all_signatures, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
es.max <- es.max_full$es.obs
es.dif <- es.dif_full$es.obs
par(mfrow=c(1,2), mar=c(4, 4, 4, 1))
plot(density(as.vector(es.max)), main="Maximum deviation from zero",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
par(mfrow=c(1,1))

# heatmap
signatures <- t(es.dif)
heatmap.plus(as.matrix(signatures), col=bluered(51), scale="none")

# add clinical annotations
clinical <- read.csv("/Users/Eric/hnscc/FromUC/HNC_Clinical_data_131104.txt", sep="\t", row.names=1)
groupings <- read.csv("/Users/Eric/hnscc/FromUC/AGI_group_allocation.txt", sep="\t", row.names=1)

rownames(signatures) <- gsub("X", "", rownames(signatures), fixed=TRUE)
rownames(groupings) <- gsub("X", "", rownames(groupings), fixed=TRUE)

# clinical$Anatomic.Site.1 has 13 categories, many only have 1-3
# groupings$group and groupings$HPV

col_mat <- c()
groupings_mat <- c()
for (i in 1:nrow(signatures))
{
  col_row <- c("","","")
  patient_id <- rownames(signatures)[i]
  clinical_row <- which(rownames(clinical)==patient_id)
  groupings_row <- which(rownames(groupings)==patient_id)
  site <- as.character(clinical$Anatomic.Site.1[clinical_row])
  group <- as.character(groupings$group[groupings_row])
  hpv <- as.character(groupings$HPV[groupings_row])
  groupings_mat <- rbind(groupings_mat, c(site, group, hpv))
  if (hpv=="pos") {col_row[3] <- "black"}
  else {col_row[3] <- "grey"}
  if (group=="BA") {col_row[2] <- "magenta"}
  else if (group=="CL") {col_row[2] <- "darkgreen"}
  else if (group=="MS") {col_row[2] <- "blue4"}
  else {col_row[2] <- "darkorange"}
  if (site=="OROPHARYNX") {col_row[1] <- "#FF0000FF"}
  else if (site=="ORAL CAVITY") {col_row[1] <- "#FFFF00FF"}
  else if (site=="SUPRAGLOTTIC") {col_row[1] <- "#00FF00FF"}
  else if (site=="HYPOPHARYNX") {col_row[1] <- "#00FFFFFF"}
  else if (site=="GLOTTIC LARYNX" || site=="SUBGLOTTIC LARYNX") {col_row[1] <- "#0000FFFF"}
  else {col_row[1] <- "#FF00FFFF"}
  col_mat <- rbind(col_mat, col_row)
}

heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none")

legend_cols <- c("#FF0000FF", "#FFFF00FF", "#00FF00FF", "#00FFFFFF", "#0000FFFF", "#FF00FFFF",
                 "magenta", "darkgreen", "blue4", "darkorange", "black", "grey")
legend_names <- c("Oropharynx", "Oral Cavity", "Supraglottic", "Hypopharynx", "Larynx",
                  "Other", "BA", "CL", "MS", "unclassified", "HPV+", "HPV-")
plot.new()
legend(0,1, legend_names, fill=legend_cols, ncol=2, cex=0.75)

### survival ###
# clusters
signatures <- scale(signatures)
d <- dist(as.matrix(signatures))
hc <- hclust(d)
mycl <- cutree(hc, h=10)

# dates
days_end <- clinical$Date.of.last.follow.up...dead.of.death
days_start <- clinical$Date.of.initial..HONC.clinic.visit...On.study.Date

# create matrix (clusters)
surv <- c()
for (i in 1:length(mycl))
{
  patient <- names(mycl[i])
  cl <- mycl[i]
  clinical_row <- which(rownames(clinical)==patient) 
  status <- as.character(clinical$Alive.Dead)[clinical_row]
  if (status=="Alive") {status <- 0}
  else if (status=="Deceased") {status <- 1}
  else (status <- "NA")
  if (nchar(as.character(days_end[clinical_row]))!=8) {years <- NA}
  else
  {
    starting_date <- as.Date(days_start[clinical_row], format='%m/%d/%y')
    ending_date <- as.Date(days_end[clinical_row], format='%m/%d/%y')
    years <- as.numeric(ending_date-starting_date)/365
  }
  surv <- rbind(surv, c(cl, years, status))
}

my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2,3), main="Survival of HNSCC by clusters", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("high", "low", "med"), fill=c(1,2,3))

# create matrix (HPV)
surv <- c()
for (i in 1:nrow(clinical))
{
  hpv <- as.character(clinical$HPV[i])
  status <- as.character(clinical$Alive.Dead)[i]
  if (status=="Alive") {status <- 0}
  else if (status=="Deceased") {status <- 1}
  else (status <- "NA")
  if (nchar(as.character(days_end[i]))!=8) {years <- NA}
  else
  {
    starting_date <- as.Date(days_start[i], format='%m/%d/%y')
    ending_date <- as.Date(days_end[i], format='%m/%d/%y')
    years <- as.numeric(ending_date-starting_date)/365
  }
  surv <- rbind(surv, c(hpv, years, status))
}

my.surv <- Surv(as.numeric(surv[,2]), as.numeric(surv[,3]))
my.fit <- survfit(my.surv ~ surv[,1])
plot(my.fit, col=c(1,2), main="Survival of HNSCC by HPV", xlab="Survival Time (years)", ylab="Survival Probability")
legend(0.5,0.4, c("neg", "pos"), fill=c(1,2))