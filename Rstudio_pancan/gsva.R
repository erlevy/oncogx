library(GSVA)

expression_dir <- "/Users/Eric/tcrseq/new/expression"

files_expression <- list.files(path=expression_dir, pattern="*.txt", full.names=T, recursive=FALSE)

expression1 <- read.csv(files_expression[1], sep="\t")[-1,]

print_lines <- function(file_path)
{
  print(summary(as.numeric(read.csv(file_path, sep="\t")[3355,-1])))
}

lapply(files_expression, print_lines)

run_gsva <- function(expression_path)
{
  cibersort <- read.csv("/Users/Eric/cibersort_signatures.txt", sep="\t")
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
  all_signature_names <- c("T-regs", "T-CD8d", "T-CD4n", "T-CD4mr", "T-CD4ma", "T-hf", "T-gd", "B-n", "B-m", "Plasma",
                           "NK-r", "NK-a", "Mono", "Macro-0", "Macro-1", "Macro-2", "Den-r", "Den-a", "Mast-r", "Mast-a",
                           "Eos", "Neut")
  names(all_signatures) <- all_signature_names
  
  data_expression <- read.csv(expression_path, sep="\t")[-1,]
  data_expression <- as.matrix(data_expression)
  rownames(data_expression) <- data_expression[,1]
  data_expression <- data_expression[,-1]
  mode(data_expression) <- 'numeric'
  gene_names <- rownames(data_expression)
  gene_names <- strsplit(gene_names, "|", fixed=TRUE)
  gene_names <- sapply(gene_names, "[", 1)
  rownames(data_expression) <- gene_names
  all_signatures <- lapply(all_signatures, as.character)
  es.dif_full <- gsva(data_expression, all_signatures, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
  es.dif <- es.dif_full$es.obs
  cancer_type <- basename(strsplit(expression_path, '[.]')[[1]][1])
  output_dir <- "/Users/Eric/tcrseq/new/gsva/"
  output_file <- paste(cancer_type, "_gsva.txt", sep="")
  gsva_out <- paste(output_dir, output_file, sep="")
  write.table(es.dif, gsva_out, sep="\t", quote=FALSE)
  print(gsva_out)
#  return(es.dif)
}

lapply(files_expression, run_gsva)
