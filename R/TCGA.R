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

row_data <- rowData(data)
colnames(row_data)
rownames(row_data)
print(row_data$gene_id)
print(row_data$gene_name)
table(row_data$gene_type)

print(metadata$definition)
length(colnames(metadata))
print(table(metadata$definition))
gene_map <- as.data.frame(row_data)[, c("gene_id", "gene_name")]
res$hgnc_symbol <- gene_map$gene_name[match(rownames(res), gene_map$gene_id)]

group <- factor(ifelse(metadata$definition %in% c("Primary solid Tumor", "Recurrent Solid Tumor"), "Tumor", "Normal"))
print(table(group))
metadata$group <- group

class(metadata$group)
str(metadata$group)


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exp_matrix,
                              colData = metadata,
                              design = ~ group)
print(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Tumor", "Normal"))
View(res)


library(ggplot2)

res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
res_df$neg_log10_padj <- -log10(res_df$padj)

res_df$significant <- ifelse(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, "Significant", "Not Significant")


ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significant)) + 
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # log2FC 기준선
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # padj 기준선
  labs(x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       title = "Volcano Plot: Tumor vs Normal",
       color = "Significance") +
  theme_minimal()

library(EnhancedVolcano)

# 상위 10개 유전자 선택
key_genes <- head(res[order(res$padj), ], 10)$log2FoldChange
names(key_genes) <- head(rownames(res[order(res$padj), ]), 10)

EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                selectLab = names(key_genes),  # 라벨링할 유전자
                title = "Volcano Plot: Tumor vs Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                drawConnectors = TRUE,         # 라벨과 점 연결선
                widthConnectors = 0.5)         # 연결선 두께



gene_name <- c()


gene_name

rownames(res) <- gene_name

key_genes <- head(res[order(res$padj), ], 10)$log2FoldChange
names(key_genes) <- head(rownames(res[order(res$padj), ]), 10)

EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "padj",
                selectLab = names(key_genes),  # 라벨링할 유전자
                title = "Volcano Plot: Tumor vs Normal",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                drawConnectors = TRUE,         # 라벨과 점 연결선
                widthConnectors = 0.5)         # 연결선 두께
