# 필요한 패키지 설치 및 로드 ===========================================
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

required_packages <- c("TCGAbiolinks", "limma", "edgeR", "tidyverse", "biomaRt", 
                       "SummarizedExperiment", "clusterProfiler", "org.Hs.eg.db", 
                       "KEGGREST", "DOSE", "DESeq2", "ggplot2", "EnhancedVolcano")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(TCGAbiolinks)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(limma)
library(edgeR)
library(tidyverse)
library(biomaRt)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(DOSE)

set.seed(123)  # 재현성 보장

# TCGA 데이터 로드 ========================================================
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

getwd()
setwd("Desktop/Study/repositories/")
getwd()

data <- GDCprepare(query, summarizedExperiment = TRUE)

exp_matrix <- assay(data)
metadata <- as.data.frame(colData(data))
row_data <- as.data.frame(rowData(data))

# 필요한 메타데이터 컬럼 선택
need_col <- c("barcode","patient","shortLetterCode", "tumor_descriptor", "sample_type", 
              "tissue_type", "age_at_diagnosis", "laterality", "treatments",
              "primary_diagnosis", "prior_treatment", "classification_of_tumor",
              "tumor_of_origin", "tobacco_smoking_status", "pack_years_smoked", 
              "gender", "race", "vital_status", "ajcc_pathologic_stage")
meta <- metadata[, need_col]

# 유전자 심볼 추가
gene_symbol <- tibble("ID" = row_data$gene_id, "Name" = row_data$gene_name)
rownames(exp_matrix) <- gene_symbol$Name[match(rownames(exp_matrix), gene_symbol$ID)]
exp_matrix <- exp_matrix[!is.na(rownames(exp_matrix)), ]

# 데이터 타입 확인
col_types <- sapply(meta, class)
print(col_types)

# 종양 단계(Stage) 정리 ===================================================
mapping <- data.frame(
  Stage = c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", 
            "Stage IIIA", "Stage IIIB", "Stage IV"),
  Numeric = c(1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.1, 3.2, 4.0)
)
meta$Stage <- mapping$Numeric[match(meta$ajcc_pathologic_stage, mapping$Stage)]
meta$Stage[is.na(meta$Stage) & meta$shortLetterCode == "NT"] <- 0.0  # NT 샘플은 0으로 설정

# TNM 정보를 기반으로 Stage 할당 (예시)
tnm_criteria <- list(
  list(T = "T1a", N = "N1", M = "M0", Stage = 2.2),
  list(T = "T2", N = "N1", M = "M0", Stage = 2.2),
  list(T = "T4", N = "N0", M = "MX", Stage = 3.1)
)

for (rule in tnm_criteria) {
  meta$Stage[is.na(meta$Stage) & meta$tissue_type == "Tumor" & 
               meta$ajcc_pathologic_t == rule$T & 
               meta$ajcc_pathologic_n == rule$N & 
               meta$ajcc_pathologic_m == rule$M] <- rule$Stage
}

# 결측값 보정
meta$Stage[is.na(meta$Stage)] <- median(meta$Stage, na.rm = TRUE)

# DEG 분석 준비 =========================================================
metadata_filtered <- meta[meta$shortLetterCode %in% c("NT", "TP"), ]
exp_matrix_filtered <- exp_matrix[, meta$shortLetterCode %in% c("NT", "TP")]


# DGEList 생성 및 정규화
group <- factor(metadata_filtered$shortLetterCode, levels = c("NT", "TP"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


dge <- DGEList(counts = exp_matrix_filtered, group = factor(metadata_filtered$shortLetterCode))
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# Voom 변환
v <- voom(dge, design = model.matrix(~ 0 + group), plot = TRUE)

# DEG 분석 수행 ==========================================================
fit <- lmFit(v, design)
contrast_matrix <- makeContrasts(Tumor_vs_Normal = TP - NT, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

# DEG 결과 확인
results <- topTable(fit2, coef = "Tumor_vs_Normal", adjust = "BH", number = Inf)
results$Significance <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Significant", "Not significant")

new_rowname <- row_data$gene_name[match(rownames(results), rownames(row_data))]
na_indices <- is.na(new_rowname)
new_rowname[na_indices] <- rownames(results)[na_indices]  # NA를 원래 이름으로 대체
new_rowname <- make.unique(new_rowname)
rownames(results) <- new_rowname
View(results)



# Volcano Plot
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("red", "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs (TP vs NT)", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

# Enhanced Volcano Plot
EnhancedVolcano(results,
                lab = results$ID,
                x = "logFC",
                y = "adj.P.Val",
                title = "Enhanced Volcano Plot (TP vs NT)",
                subtitle = "TCGA-LUAD DEG Analysis",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 4.0,
                col = c("black", "blue", "green", "red"),
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 5)

immune_markers <- read.csv("LM22.txt", sep = "\t")
immune_markers <- immune_markers$Gene.symbol
print(immune_markers)



# 상위 DEG Boxplot ====================================================
top_genes <- head(results[order(results$adj.P.Val), ], 10)$ID 
exp_long <- as.data.frame(t(exp_matrix[top_genes, ]))
exp_long$Sample <- rownames(exp_long)
exp_long <- reshape2::melt(exp_long, id.vars = "Sample")
colnames(exp_long) <- c("Sample", "Gene", "Expression")
exp_long$Group <- metadata_filtered$shortLetterCode[match(exp_long$Sample, rownames(metadata_filtered))]
exp_long$Expression <- log2(exp_long$Expression + 1)

Boxplot_gene <- function(df, g) {
  ggplot(subset(df, Gene == g), aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(0.2), alpha = 0.5) +
    theme_minimal() +
    scale_fill_manual(values = c("NT" = "blue", "TP" = "red")) +
    labs(title = paste("Expression of", g, "(NT vs TP)"), y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
# TME ==============================================================


# UMAP =============================================================
if (!requireNamespace("uwot", quietly = TRUE)) 
  install.packages("uwot")

library(uwot)
umap_results <- umap(t(exp_matrix_filtered), 
                     n_neighbors = 15, 
                     min_dist = 0.1, 
                     metric = "euclidean")
umap_df <- data.frame(
  UMAP1 = umap_results[,1],
  UMAP2 = umap_results[,2],
  Group = metadata_filtered$shortLetterCode
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("NT" = "blue", "TP" = "red")) +
  labs(title = "UMAP Clustering of TCGA-LUAD Samples", x = "UMAP1", y = "UMAP2")

# Survival =========================================================





