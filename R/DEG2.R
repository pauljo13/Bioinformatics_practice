# DEG 2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("DEGseq")
BiocManager::install("DESeq2")

install.packages("pheatmap")

library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)

gse <- getGEO("GSE15824", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse <- gse[[1]]

exp <- exprs(gse)
View(exp)

if (max(exp) > 50) {
  exp <- log2(exp + 1)
}
View(exp)

pheno <- pData(gse)
View(pheno)
pheno$condition <- ifelse(grepl("normal", pheno$title, ignore.case = TRUE), 
                               "Normal", 
                               "GBM")
View(pheno)

colData <- data.frame(condition = pheno$condition, row.names = pheno$geo_accession)
colData
colData$condition <- factor(colData$condition, levels = c("Normal", "GBM"))

exp <- exp[, match(rownames(colData), colnames(exp))]

all(colnames(exp) == rownames(colData))

dds <- DESeqDataSetFromMatrix(countData = round(exp),
                              colData = colData,
                              design = ~ condition)

