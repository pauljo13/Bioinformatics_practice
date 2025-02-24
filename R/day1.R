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
pheno <- attr(attr(ges[[1]], "phenoData"), "data")
exprs <- attr(ges[[1]], "assayData")[["exprs"]]

dim(pheno)
colnames(pheno)
dim(exprs)
colnames(exprs)

colnames(exprs) <- pheno[colnames(exprs),"title"]
colnames(exprs)

rownames(exprs)
nan_count <- sum(is.na(rownames(exprs)) | rownames(exprs) == "NaN")
print(nan_count)

setwd("/Users/knu_cgl1/Desktop/Study/repositories/GSE15824")
write.table(exprs, file = "GSE15824_expression_original.txt", sep = "\t")

ges[[1]]
gpl <- getGEO("GPL570")
gpl_table <- Table(gpl)
colnames(gpl_table)
gpl_table["Gene Symbol"]
gpl_table["ID"]


original_rownames <- rownames(exprs)
length(original_rownames)

matched_indices <- match(original_rownames, gpl_table$ID)
length(matched_indices)
gene_symbols <- gpl_table$`Gene Symbol`[matched_indices]
length((gene_symbols))

# 매칭되지 않은 값들 확인
unmatched_indices <- is.na(gene_symbols)  # NA인 값 찾기
unmatched_genes <- original_rownames[unmatched_indices]  # 매칭되지 않은 원래 rownames 저장

# NA가 있는 경우, 원래의 rownames 사용
gene_symbols[unmatched_indices] <- original_rownames[unmatched_indices]

# rownames 변경
rownames(exprs) <- gene_symbols



getwd()
setwd("/Users/knu_cgl1/Desktop/Study/repositories/GSE15824")
write.table(exprs, file = "GSE15824_expression.txt", sep = "\t")
write.table(pheno, file = "GSE15824_phenotype.txt", sep = "\t")
write.table(gpl_table, file = "GPL570.txt", sep = "\t")
