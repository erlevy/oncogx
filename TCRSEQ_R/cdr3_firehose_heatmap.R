library(gplots)

load("cibersort_df.Rdata")
cdr3 <- read.csv("cdr3_results_with_gsva.csv")
cdr3$cdr3_grouping <- gsub("zero", "low", cdr3$cdr3_grouping, fixed=TRUE)
cdr3 <- cdr3[complete.cases(cdr3),]

all_meds <- matrix(nrow=6, ncol=22)
for (i in 1:ncol(all_meds))
{
  signature_start_index <- 19
  cdr3$signature_all <- cdr3[,19+i]
  
  split_high_low <- split(cdr3, interaction(cdr3$clinical_group, cdr3$cdr3_group))
  hr_positive_high <- split_high_low[[1]]
  hr_positive_low <- split_high_low[[5]]
  her2_positive_high <- split_high_low[[2]]
  her2_positive_low <- split_high_low[[6]]
  tnbc_high <- split_high_low[[4]]
  tnbc_low <- split_high_low[[8]]

  # signature_all
  her2_wilcox_data_signature_all <- rbind(cbind(her2_positive_high$signature_all, 1), cbind(her2_positive_low$signature_all, 0))
  hr_wilcox_data_signature_all <- rbind(cbind(hr_positive_high$signature_all, 1), cbind(hr_positive_low$signature_all, 0))
  tnbc_wilcox_data_signature_all <- rbind(cbind(tnbc_high$signature_all, 1), cbind(tnbc_low$signature_all, 0))
  
  # with normalization
#  her2_high_med <- median(her2_positive_high$signature_all)/median(cdr3$signature_all)
#  her2_low_med <- median(her2_positive_low$signature_all)/median(cdr3$signature_all)
#  hr_high_med <- median(hr_positive_high$signature_all)/median(cdr3$signature_all)
#  hr_low_med <- median(hr_positive_low$signature_all)/median(cdr3$signature_all)
#  tnbc_high_med <- median(tnbc_high$signature_all)/median(cdr3$signature_all)
#  tnbc_low_med <- median(tnbc_low$signature_all)/median(cdr3$signature_all)
  
  #without normalization
  her2_high_med <- median(her2_positive_high$signature_all)
  her2_low_med <- median(her2_positive_low$signature_all)
  hr_high_med <- median(hr_positive_high$signature_all)
  hr_low_med <- median(hr_positive_low$signature_all)
  tnbc_high_med <- median(tnbc_high$signature_all)
  tnbc_low_med <- median(tnbc_low$signature_all)
  
  meds <- c(her2_high_med, her2_low_med, hr_high_med, hr_low_med, tnbc_high_med, tnbc_low_med)
  all_meds[,i] <- meds
}


all_signature_names <- c("T regs", "T-CD8", "T-CD4n", "T-CD4mr", "T-CD4ma", "T-hf", "T-gd", "B-n", "B-m", "Plasma",
                         "NK-r", "NK-a", "Mono", "Macro-0", "Macro-1", "Macro-2", "Den-r", "Den-a", "Mast-r", "Mast-a",
                         "Eos", "Neut")
all_row_names <- c("HER2+ high", "HER2+ low", "HR+ high", "HR+ low", "TNBC high", "TNBC low")

colnames(all_meds) <- all_signature_names
rownames(all_meds) <- all_row_names

heatmap.2(all_meds, trace="none", dendrogram='none', Rowv='none', Colv='none', col=greenred(16), margins = c(6,9))


# trying out heatmaps
heatmap.2()
