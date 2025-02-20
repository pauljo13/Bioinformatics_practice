# Study Bioinformatics Day 1
# install package and load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DEGseq")
BiocManager::install("GEOquery")

library(GEOquery)
library(tibble)
library(DEGseq)
library(ggplot2)
library(dplyr)
library(splines)

# data load
ges <- getGEO("GSE15824")
ges[[1]]
type(ges[[1]])
attr(ges[[1]], "phenoData")
pheno <- attr(attr(ges[[1]], "phenoData"), "data")
exprs <- attr(ges[[1]], "assayData")[["exprs"]]

dim(pheno)
colnames(pheno)
dim(exprs)
colnames(exprs)

colnames(exprs) <- pheno[colnames(exprs),"title"]
colnames(exprs)

ges[[1]]
gpl <- getGEO("GPL570")
gpl_table <- Table(gpl)
colnames(gpl_table)
gpl_table["Gene Symbol"]
gpl_table["ID"]

gene_symbols <- gpl_table$`Gene Symbol`[match(rownames(exprs), gpl_table$ID)]
colnames(exprs)
rownames(exprs) <- gene_symbols
rownames(exprs)

getwd()
setwd("/Users/knu_cgl1/Desktop/Study/repositories/GSE15824")
write.table(exprs, file = "GSE15824_expression.txt", sep = "\t")
write.table(pheno, file = "GSE15824_phenotype.txt", sep = "\t")
