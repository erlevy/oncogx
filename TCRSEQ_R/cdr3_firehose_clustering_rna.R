library(OIsurv)
library(cgdsr)
library(cluster)
library(fpc)
library(heatmap.plus)

### messing around with lymphocytes
results <- read.csv("cdr3_results_rna_with_gsva.csv")
results <- na.omit(results)
# get rid of her2+
results <- results[which(results$clinical_group!="her2+"),]
# get rid of "other"
results <- results[which(results$clinical_group!="other"),]

# isole gene signature values
signatures <- results[,20:41]
rownames(signatures) <- results[,1]

# normalize?
signatures <- scale(signatures)

# Determine number of clusters
wss <- (nrow(signatures)-1)*sum(apply(signatures,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(signatures, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# K-Means Cluster Analysis
fit <- kmeans(signatures, 4)
# get cluster means 
aggregate(signatures,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(results, fit$cluster)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
clusplot(signatures, fit$cluster, color=TRUE, shade=TRUE, labels=0, lines=0)
# Centroid Plot against 1st 2 discriminant functions
plotcluster(signatures, fit$cluster)

# hierarchical cluster analysis
d <- dist(as.matrix(signatures))
hc <- hclust(d)
plot(hc, labels=FALSE)
mycl <- cutree(hc, h=10)
#mycl <- cutree(hc, h=max(hc$height/1.5))

results_hc <- cbind(results, mycl)

col_mat <- matrix(nrow=nrow(results), ncol=3)
exome_freq <- results$exome_imseq/results$exome_reads
rna_freq <- results$rna_imseq/results$rna_reads
exome_med <- median(exome_freq[which(exome_freq!=0)])
rna_med <- median(rna_freq[which(rna_freq!=0)])
for (i in 1:nrow(col_mat))
{
  group <- results$clinical_group[i]
  exome <- exome_freq[i]
  rna <- rna_freq[i]
  if (exome<=exome_med) {col_mat[i,2]="white"}
  else {col_mat[i,2]="black"}
  if (rna<=rna_med) {col_mat[i,3]="white"}
  else {col_mat[i,3]="purple"}
  if (group=="her2-") {col_mat[i,1]="red"}
  else {col_mat[i,1]="green"}
}

heatmap.plus(as.matrix(signatures), col=bluered(51), RowSideColors=col_mat, scale="none")
