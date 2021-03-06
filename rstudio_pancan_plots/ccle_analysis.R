library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(gplots)
library(heatmap.plus)
library(reshape)
library(dplyr)
library(ggplot2)
library(plotrix)

results <- read.csv("/Users/Eric/tcrseq/new/processed/ccle_imseq_results_7-22-2016.txt", sep="\t")
clonotypes <- read.csv("/Users/Eric/tcrseq/new/processed/ccle_clonotypes.txt", sep="\t", header=TRUE)

results_pos <- filter(results, imseq > 0)
