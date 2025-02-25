# DEG
# library insall 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

library(limma)

getwd()
setwd("/Users/bumsoojoe/Documents/Obsidian Vault/programming/DATA")
getwd()

df <- read.csv("df_gene.csv", row.names = 1)
colnames(df)

gsms <- 