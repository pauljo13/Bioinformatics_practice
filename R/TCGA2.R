# DEG used TCGA data
# istall packages and load library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManger")

BiocManager::install("TCGAbiolinks")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("tidyverse")
BiocManager::install("biomaRt")


library(TCGAbiolinks)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

set.seed(123)  # 재현성 보장
library(limma)
library(edgeR)
library(tidyverse)
library(biomaRt)

# data load
# TCGA-LUAD 데이터를 검색
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM"
)

data <- GDCprepare(
  query,
  summarizedExperiment = TRUE
)

exp_matrix <- assay(data)
metadata <- colData(data)
row_data <- rowData(data)

View(exp_matrix)
View(metadata)
View(row_data)

gene_symbol <- tibble("ID" = row_data$gene_id,"Name" = row_data$gene_name,"Type" = row_data$gene_type)
gene_symbol


new_rowname <- gene_symbol$Name[match(rownames(exp_matrix), gene_symbol$ID)]
rownames(exp_matrix) <- new_rowname
exp_matrix <- exp_matrix[!is.na(rownames(exp_matrix)), ]

dim(exp_matrix)

# DEG analysis
dim(exp_matrix)
table(metadata$shortLetterCode)

metadata_filtered <- metadata[metadata$shortLetterCode %in% c("NT", "TP"), ]
exp_matrix_filtered <- exp_matrix[, metadata$shortLetterCode %in% c("NT", "TP")]
table(metadata_filtered$shortLetterCode)

# DGElist 
dge <- DGEList(counts = exp_matrix_filtered)

# TMM normalization
dge <- calcNormFactors(dge)

group <- factor(metadata_filtered$shortLetterCode, levels = c("NT", "TP"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Voom
v <- voom(dge, design = design, plot = TRUE)

# model fit
fit <- lmFit(v, design)

# Contrasts
contrast_matrix <- makeContrasts(Tumor_vs_Normal = TP - NT, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

View(fit2)

topTable(fit2, coef = "Tumor_vs_Normal", adjust = "BH", number = 20)

# Volcano Plot
deg_results$Significance <- ifelse(
  deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1, "Significant", "Not significant"
  )

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("red", "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes (TP vs NT)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

EnhancedVolcano(deg_results,
                lab = rownames(deg_results),           # 유전자 이름
                x = "logFC",                           # x축: Log Fold Change
                y = "adj.P.Val",                       # y축: Adjusted P-value
                xlab = bquote(~Log[2]~" Fold Change"), # x축 라벨
                ylab = bquote(~-Log[10]~" Adjusted P-value"), # y축 라벨
                title = "Enhanced Volcano Plot (TP vs NT)", # 제목
                subtitle = "TCGA-LUAD DEG Analysis",  # 부제목
                pCutoff = 0.05,                        # p-value 기준
                FCcutoff = 1,                          # Log2 FC 기준 (1 이상이면 의미 있음)
                pointSize = 2.0,                        # 점 크기
                labSize = 4.0,                          # 라벨 크기
                colAlpha = 0.75,                        # 점 투명도
                legendPosition = "right",               # 범례 위치
                legendLabSize = 12,                     # 범례 폰트 크기
                legendIconSize = 4.0,                    # 범례 아이콘 크기
                col = c("black", "blue", "green", "red"), # 점 색상 설정
                drawConnectors = TRUE,                  # 유전자명과 점을 선으로 연결
                widthConnectors = 0.5                    # 선 굵기
)


head(rownames(deg_results))
head(rownames(exp_matrix_filtered))

all(rownames(exp_matrix_filtered) %in% rownames(deg_results))
all(rownames(deg_results) %in% rownames(exp_matrix_filtered))

head(rownames(v))
