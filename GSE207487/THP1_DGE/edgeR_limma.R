library(edgeR)
library(limma)
library(Glimma)
library(tidyverse)
library(Homo.sapiens)
library(RColorBrewer)

# Data Organization ----------------------------
counts_matrix <- read.csv("../data/raw_counts_GRCh38_GREIS.csv") %>% 
  dplyr::select(c(1:5, 7:9)) %>% # compare uninfected to WT infected cells, triplicate
  column_to_rownames(var = "X")

group <- as.factor(c(rep("UN", 3), rep("WT", 3)))
x <- DGEList(counts_matrix, group = group)
x$samples

# first column contains old gene symbols, updated version for HGNC ids
geneid <- rownames(x) # Entrez gene ids from Salmon
genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM"), 
                keytype = "ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),]
head(genes)

x$genes <- genes
x

# Normalization and Filtering --------------------------------
cpm <- cpm(x)
lcpm <- cpm(x, log=T) # log2(CPM + 2/L) to add "prior" value
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6

table(rowSums(x$counts==0)==6) # zero counts across all samples

# Figure 1 - The Density of log-CPM values before/after filtering for low genes
samplenames <- colnames(x$counts)
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")

# A - Raw Data
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# Filter
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# B - Filtered Data
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# calculate normalization with TMM (Oshlack 2010)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# Figure 2 MDS Plot
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
plotMDS(lcpm, labels=group)
title(main="Sample groups MDS")


