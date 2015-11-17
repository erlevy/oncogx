library(OIsurv)
library(cgdsr)
library(GSVA)

expression <- read.csv("tcga/expression/raw/Firehose/BRCA_normalized.csv")
#expression <- read.csv("BRCA_expression_firehose_patient_ids.csv")
cdr3 <- read.csv("cdr3_results_with_groups.txt", sep="\t")

# tutorial
p <- 20000 ## number of genes
n <- 30 ## number of samples
nGS <- 100 ## number of gene sets
min.sz <- 10 ## minimum gene set size
max.sz <- 100 ## maximum gene set size
X <- matrix(rnorm(p*n), nrow=p, dimnames=list(1:p, 1:n))
dim(X)
gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE)) ## sample gene set sizes
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets
es.max <- gsva(X, gs, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)$es.obs
es.dif <- gsva(X, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)$es.obs
par(mfrow=c(1,2), mar=c(4, 4, 4, 1))
plot(density(as.vector(es.max)), main="Maximum deviation from zero",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
plot(density(as.vector(es.dif)), main="Difference between largest\npositive and negative deviations",
     xlab="GSVA score", lwd=2, las=1, xaxt="n", xlim=c(-0.75, 0.75), cex.axis=0.8)
axis(1, at=seq(-0.75, 0.75, by=0.25), labels=seq(-0.75, 0.75, by=0.25), cex.axis=0.8)
par(mfrow=c(1,1))

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
expression <- as.matrix(expression)
rownames(expression) <- expression[,1]
expression <- expression[,-1]
mode(expression) <- 'numeric'
gene_names <- rownames(expression)
gene_names <- strsplit(gene_names, "|", fixed=TRUE)
gene_names <- sapply(gene_names, "[", 1)
rownames(expression) <- gene_names
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

# use uuids instead of barcodes (uuid_to_barcode.R)
#colnames(es.max) <- uuids
#colnames(es.dif) <- uuids

# add gsva results to cdr3 results table
#es.dif_col <- c()
#es.dif_all <- gsub(".", "-", colnames(es.dif), fixed=TRUE)
#es.dif_all <- gsub("X", "", es.dif_all, fixed=TRUE)
#for (i in 1:nrow(cdr3))
#{
#  rna <- as.character(cdr3$patient_id[i])
#  es.dif_index <- which(rna==es.dif_all)
#  if (length(es.dif_index)>0)
#  {
#    es.dif_value <- es.dif[,es.dif_index]
#    es.dif_col <- rbind(es.dif_col, es.dif_value)
#  }
#  else
#  {
#    es.dif_col <- rbind(es.dif_col, rep(NA, 22))
#  }
#}
#rownames(es.dif_col) <- seq(1,nrow(cdr3))
#cdr3_new <- cbind(cdr3, es.dif_col)
#write.csv(cdr3_new, "cdr3_results_with_gsva.csv", quote=FALSE, row.names=FALSE)