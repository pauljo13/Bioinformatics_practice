# DEG
# library insall 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("DEGseq")
BiocManager::install("GEOquery")

library(GEOquery)
library(tibble)
library(DEGseq)
library(ggplot2)
library(dplyr)
library(splines)
library(limma)

getwd()
setwd("/Users/bumsoojoe/Documents/Obsidian Vault/programming/DATA")
getwd()

df <- read.csv("df_gene.csv")
colnames(df)

set.seed(123)
gset <- getGEO("GSE15824", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

# 1 : GBM, 0 : Normal
sml <- c()
for (i in colnames(df)) {
  if (grepl("normal", i)) {
    sml <- c(sml, "0")
  } else if (grepl("cell", i)) {
    sml <- c(sml, "X")
  } else {
    sml <- c(sml, "1")
  }
}; sml
sel <- which(sml != "X"); sel
sml <- sml[sel]
gset <- gset[,sel]

ex <- as.matrix(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

gs <- factor(sml); gs
groups <- make.names(c("norml", "GBM"))
levels(gs) <- groups; gs
length(gs)
length(colnames(gset))
gset$group <- gs
gset$group
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)

cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)
deg_results <- topTable(fit2, adjust = "fdr", sort.by = "B", number = lnf)
