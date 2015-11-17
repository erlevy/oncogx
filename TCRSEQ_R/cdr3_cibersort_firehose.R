library(OIsurv)
library(cgdsr)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Get available case lists (collection of samples) for brca tcga
tcga_brca = getCancerStudies(mycgds)[16,1]
brca_caselist = getCaseLists(mycgds,tcga_brca)[8,1] # tumors with mRNA data (1100)
brca_clinical = getClinicalData(mycgds,brca_caselist)
brca_genetic = getGeneticProfiles(mycgds,tcga_brca)[4,1] # RNA-seq Z-scores

# random gene ids
ids <- read.csv("tcga_gene_ids.txt", sep="\t")
ids <- as.matrix(ids)[,1]
ids_sample <- sample(ids, 50)

# cibersort gene signatures
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
all_row_names <- c("her2+ high/low", "hr+ high/low", "tnbc high/low",
                   "her2+/hr+ high", "her2+/tnbc high", "hr+/tnbc high",
                   "her2+/hr+ all", "her2+/tnbc all", "hr+/tnbc all")
#control <- ids_sample
# cdr3 results
results <- read.csv("cdr3_results_with_groups.txt", sep="\t")

cdr3_gene_df <- function(genes, results)
{
  brca_genes = getProfileData(mycgds, genes, brca_genetic, brca_caselist)
  if (sum(is.na(brca_genes[1,]))>0) {brca_genes <- brca_genes[,-which(is.na(brca_genes[1,]))]}
  gene_signature <- rowSums(brca_genes)
  cdr3_freq <- results$exome_imseq/results$exome_reads
  cor_table <- c()
  cdr3_all <- c()
  signature_all <- c()
  groups_all <- c()
  sample_genes <- c()
  genes_all <- c()
  cdr3_group <- c()
  cdr3_med <- median(cdr3_freq[which(cdr3_freq!=0)])
  for (i in 1:nrow(results))
  {
    sample_id <- as.character(results$sample_id[i])
    sample_id <- substr(sample_id, 1, 15)
    sample_id <- gsub("-", ".", sample_id)
    sample_signature <- gene_signature[sample_id]
    sample_genes <- brca_genes[sample_id,]
    sample_cdr3 <- as.numeric(results$exome_imseq[i])/as.numeric(results$exome_reads[i])
    sample_group <- as.character(results$clinical_group[i])
    if (is.na(sample_signature)) {next}
    else
    {
      cor_table <- rbind(cor_table, c(sample_cdr3, as.numeric(sample_signature)))
      cdr3_all <- c(cdr3_all, sample_cdr3)
      signature_all <- c(signature_all, sample_signature)
      groups_all <- c(groups_all, sample_group)
      genes_all <- rbind(genes_all, sample_genes)
      if (sample_cdr3 < cdr3_med) {cdr3_group <- c(cdr3_group, "low")}
      else {cdr3_group <- c(cdr3_group, "high")}
    }
  }
  df <- data.frame(cdr3_all, genes_all, groups_all, cdr3_group, signature_all)
  return(df)
}

all_pvals <- matrix(nrow=9, ncol=22)
all_dfs <- list()

for (i in 1:length(all_signatures))
{
  genes <- as.matrix(all_signatures[[i]])[,1]
  df <- cdr3_gene_df(genes, results)
  all_dfs[[i]] <- df
}

for (i in 1:length(all_signatures))
{
#  genes <- as.matrix(all_signatures[[i]])[,1]
  
#  df <- cdr3_gene_df(genes, results)
 
  df <- all_dfs[[i]]
  
  split_high_low <- split(df, interaction(df$groups_all, df$cdr3_group))
  hr_positive_high <- split_high_low[[1]]
  hr_positive_low <- split_high_low[[5]]
  her2_positive_high <- split_high_low[[2]]
  her2_positive_low <- split_high_low[[6]]
  tnbc_high <- split_high_low[[4]]
  tnbc_low <- split_high_low[[8]]
  
  graph_title <- paste(all_signature_names[i], "expression by CDR3 high/low groups")
  boxplot(her2_positive_high$signature_all, her2_positive_low$signature_all, hr_positive_high$signature_all, hr_positive_low$signature_all, tnbc_high$signature_all, tnbc_low$signature_all, 
          names=c("HER2+ high", "HER2+ low", "HR+ high", "HR+ low", "TNBC high", "TNBC low"), ylab="Expression Z-score", main=graph_title)
  
  # signature_all
  her2_wilcox_data_signature_all <- rbind(cbind(her2_positive_high$signature_all, 1), cbind(her2_positive_low$signature_all, 0))
  hr_wilcox_data_signature_all <- rbind(cbind(hr_positive_high$signature_all, 1), cbind(hr_positive_low$signature_all, 0))
  tnbc_wilcox_data_signature_all <- rbind(cbind(tnbc_high$signature_all, 1), cbind(tnbc_low$signature_all, 0))
  
  wilcox_pvals <- rep(1,9)
  wilcox_pvals[1] <- wilcox.test(her2_wilcox_data_signature_all[,1]~her2_wilcox_data_signature_all[,2])$p.val
  wilcox_pvals[2] <- wilcox.test(hr_wilcox_data_signature_all[,1]~hr_wilcox_data_signature_all[,2])$p.val
  wilcox_pvals[3] <- wilcox.test(tnbc_wilcox_data_signature_all[,1]~tnbc_wilcox_data_signature_all[,2])$p.val
  
  wilcox_pvals[4] <-  wilcox.test(her2_positive_high$signature_all, hr_positive_high$signature_all)$p.val
  wilcox_pvals[5] <-  wilcox.test(her2_positive_high$signature_all, tnbc_high$signature_all)$p.val
  wilcox_pvals[6] <-  wilcox.test(hr_positive_high$signature_all, tnbc_high$signature_all)$p.val
  
  wilcox_pvals[7] <-  wilcox.test(her2_wilcox_data_signature_all[,1], hr_wilcox_data_signature_all[,1])$p.val
  wilcox_pvals[8] <-  wilcox.test(her2_wilcox_data_signature_all[,1], tnbc_wilcox_data_signature_all[,1])$p.val
  wilcox_pvals[9] <-  wilcox.test(hr_wilcox_data_signature_all[,1], tnbc_wilcox_data_signature_all[,1])$p.val
  
  all_pvals[,i] <- wilcox_pvals
}

colnames(all_pvals) <- all_signature_names
rownames(all_pvals) <- all_row_names
write.csv(all_pvals, "all_cibersort_pvals.csv", quote=FALSE)

# check with ucsc expression data
ucsc <- read.table("BRCA_expression_all.txt", sep="\t", )
load("cibersort_df.Rdata")
cibersort <- read.csv("cibersort_signatures.txt", sep="\t")
tregs <- cibersort$genes[which(cibersort$T.cells.regulatory..Tregs.==1)]
treg_expression <- ucsc[,tregs]
treg_cbioportal <- all_dfs[[1]]
cbioportal_patients <- rownames(treg_cbioportal)

sample_genes <- c("CD2", "CD28", "CD4", "FOXP3", "GZMM", "NTN3")
treg_expression <- ucsc[,sample_genes[4], drop=FALSE]
signature_pairs <- c()
for (i in 1:nrow(treg_expression))
{
  ucsc_row <- treg_expression[i,,drop=FALSE]
  id <- rownames(ucsc_row)
  id_match <- which(cbioportal_patients == id)
  if (length(id_match ==1))
  {
    cbioportal_signature <- treg_cbioportal[id_match,sample_genes[1]]
    ucsc_signature <- sum(ucsc_row)
    signature_pairs <- rbind(signature_pairs, c(ucsc_signature, cbioportal_signature))
  }
}

plot(signature_pairs)
abline(lm(signature_pairs[,2]~signature_pairs[,1]))
cor.test(signature_pairs[,1], signature_pairs[,2])
