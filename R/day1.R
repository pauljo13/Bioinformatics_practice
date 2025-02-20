# Study Bioinformatics Day 1
# install package and load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DEGseq")
BiocManager::install("GEOquery")

library(GEOquery)

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
dim(exprs)


