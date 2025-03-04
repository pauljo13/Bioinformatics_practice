# User TCGA data
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManger")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
); query

getwd()
setwd("/Users/knu_cgl1/Desktop/Study/repositories")
getwd()

GDCdownload(
  query = query,
  method = "api",
  directory = "GDCdata"
)

data <- GDCprepare(
  query,
  summarizedExperiment = TRUE
)

head(assay(data))
colData(data)
rowData(data)

exp_matrix <- assay(data)
head(exp_matrix[,1:5])
View(exp_matrix)

metadata <- colData(data)
View(metadata)

print(metadata$definition)
length(colnames(metadata))
print(table(metadata$definition))

group <- ifelse(metadata$definition %in% c("Primary solid Tumor", "Recurrent Solid Tumor"), "Tumor", "Normal")
print(table(group))
metadata$group <- group

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = metadata,
                              design = ~ group)
