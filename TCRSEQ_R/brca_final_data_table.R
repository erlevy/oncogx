library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)

### January 12, 2016 ###
# Purpose is to generate final tables and figures for the BRCA TCRSEQ paper

# Input: final results table to be used in paper
results1 <- read.csv("cdr3_results.txt", sep="\t")
results2 <- read.csv("cdr3_results_with_groups.txt", sep="\t")
results3 <- 