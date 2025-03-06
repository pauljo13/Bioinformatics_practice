# DEG used TCGA data
# install packages and load library ===========================================
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

# data load ===================================================================
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


# EDA ==========================================================================
View(exp_matrix)
View(metadata)
View(row_data)

dim(exp_matrix)
dim(row_data)

meta <-as.data.frame(metadata)
need_col <- c("barcode","patient","shortLetterCode", "tumor_descriptor", "sample_type", 
              "tissue_type", "age_at_diagnosis", "laterality", "treatments",
              "tissue_or_organ_of_origin", "primary_diagnosis", "prior_treatment", "classification_of_tumor",
              "icd_10_code", "tumor_of_origin", "sites_of_involvement", "tobacco_smoking_quit_year",
              "tobacco_smoking_status", "pack_years_smoked", "gender", "race", "vital_status",
              "age_at_index", "days_to_death", "primary_site", "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m")
meta <- meta[,need_col]
dim(meta)
View(meta)


gene_symbol <- tibble("ID" = row_data$gene_id,"Name" = row_data$gene_name,"Type" = row_data$gene_type)
gene_symbol

new_rowname <- gene_symbol$Name[match(rownames(exp_matrix), gene_symbol$ID)]
rownames(exp_matrix) <- new_rowname
exp_matrix <- exp_matrix[!is.na(rownames(exp_matrix)), ]

# data type
col_types <- sapply(meta, class)
print(col_types)
int <- names(col_types[col_types == "integer"])
chr <- names(col_types[col_types == "character"])
lis <- names(col_types[col_types == "list"])
numer <- names(col_types[col_types == "numeric"]) 

### data sample type
freq <- table(meta$shortLetterCode)
percent <- round(prop.table(freq) * 100, 2)  # prop.table로 비율 계산 후 퍼센트로 변환
labels_with_percent <- paste(names(freq), " (", percent, "%)", sep = "")
pie(freq, labels = labels_with_percent, col = c("red", "blue", "green"), main="shortLetterCode")

## Tumor stage
stage <- table(meta$ajcc_pathologic_stage)
print(stage)
print(names(stage))

print(sum(is.na(meta$ajcc_pathologic_stage)))
print(round((sum(is.na(meta$ajcc_pathologic_stage)) / dim(meta)[1]) * 100, 2))
mapping <- data.frame(
  Stage = c("Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", 
            "Stage IIIA", "Stage IIIB", "Stage IV"),
  Numeric = c(1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.1, 3.2, 4.0)
)
meta$Stage <- mapping$Numeric[match(meta$ajcc_pathologic_stage, mapping$Stage)]

sum(is.na(meta$Stage))
table(meta$shortLetterCode)
mean(meta[meta$shortLetterCode == "TP", ]$Stage, na.rm = TRUE)

### Stage 결측치 처리
stNa <- meta[is.na(meta$Stage),]
tnm_not_na <- stNa[!is.na(stNa$ajcc_pathologic_t) & 
                     !is.na(stNa$ajcc_pathologic_n) & 
                     !is.na(stNa$ajcc_pathologic_m), ]
tnm_not_na[,c("Stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", 
              "tissue_type", "shortLetterCode")]

### Stage NA 중 NT인 값 처리
dim(stNa[stNa$shortLetterCode == "NT",])
sum(is.na(meta$Stage))
meta$Stage[meta$shortLetterCode == "NT" & is.na(meta$ajcc_pathologic_stage)] <- 0.0
sum(is.na(meta$Stage))

### NA인 Tumor 샘플에 TNM 기반으로 Numeric 채우기
sum(is.na(meta$Stage))
meta$Stage[is.na(meta$Stage) & meta$tissue_type == "Tumor" & 
               !is.na(meta$ajcc_pathologic_t) & !is.na(meta$ajcc_pathologic_n) & 
               !is.na(meta$ajcc_pathologic_m) & 
               meta$ajcc_pathologic_t == "T1a" & 
               meta$ajcc_pathologic_n == "N1" & meta$ajcc_pathologic_m == "M0"] <- 2.2  # Stage IIB
meta$Stage[is.na(meta$Stage) & meta$tissue_type == "Tumor" & 
               !is.na(meta$ajcc_pathologic_t) & !is.na(meta$ajcc_pathologic_n) & 
               !is.na(meta$ajcc_pathologic_m) & 
               meta$ajcc_pathologic_t == "T2" & 
               meta$ajcc_pathologic_n == "N1" & meta$ajcc_pathologic_m == "M0"] <- 2.2  # Stage IIB
meta$Stage[is.na(meta$Stage) & meta$tissue_type == "Tumor" & 
               !is.na(meta$ajcc_pathologic_t) & !is.na(meta$ajcc_pathologic_n) & 
               !is.na(meta$ajcc_pathologic_m) & 
               meta$ajcc_pathologic_t == "T4" & 
               meta$ajcc_pathologic_n == "N0" & meta$ajcc_pathologic_m == "MX"] <- 3.1  # Stage IIIA (M0 가정)
sum(is.na(meta$Stage))

tnm_not_na <- stNa[!is.na(stNa$ajcc_pathologic_t) | 
                     !is.na(stNa$ajcc_pathologic_n) | 
                     !is.na(stNa$ajcc_pathologic_m), ]
tnm_not_na[,c("Stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m", 
              "tissue_type", "shortLetterCode")]

meta$Stage[is.na(meta$Stage) & meta$tissue_type == "Tumor" & 
               meta$ajcc_pathologic_t == "T1a" & 
               meta$ajcc_pathologic_n == "NX" & is.na(meta$ajcc_pathologic_m)] <- 1.1  # Stage IA (NX = N0 가정)
meta$Stage[is.na(meta$Stage) & meta$tissue_type == "Tumor" & 
               meta$ajcc_pathologic_t == "T2b" & 
               is.na(meta$ajcc_pathologic_n) & is.na(meta$ajcc_pathologic_m)] <- 2.1  # Stage IIA (N0, M0 가정)
sum(is.na(meta$Stage)) / dim(meta)[1]
stNa <- meta[is.na(meta$Stage),]

tumor_median <- median(meta$Stage[meta$tissue_type == "Tumor"], na.rm = TRUE)
cat("Tumor 샘플의 중앙값:", tumor_median, "\n")

meta$Stage[is.na(meta$Stage)] <- tumor_median
sum(is.na(meta$Stage))

### Stage 이상치 확인
NT_stage <- meta[meta$shortLetterCode == "NT" & meta$Stage != 0.0,]
dim(NT_stage)
NT_all_tnm_na <- NT_stage[is.na(NT_stage$ajcc_pathologic_t) & 
                     is.na(NT_stage$ajcc_pathologic_n) & 
                     is.na(NT_stage$ajcc_pathologic_m), ]
NT_tnm_na <- NT_stage[is.na(NT_stage$ajcc_pathologic_t) | 
                            is.na(NT_stage$ajcc_pathologic_n) | 
                            is.na(NT_stage$ajcc_pathologic_m), ]
dim(NT_all_tnm_na)
dim(NT_tnm_na)

#### NT_stage와 동일 환자의 TP 샘플 추출
NT_patients <- unique(NT_stage$patient)
TP_stage <- meta[meta$shortLetterCode == "TP" & meta$patient %in% NT_patients, ]

#### 환자별 Stage 비교
library(dplyr)
stage_comparison <- meta %>%
  filter(patient %in% NT_patients) %>%
  group_by(patient) %>%
  summarise(
    NT_Stage = unique(Stage[shortLetterCode == "NT"]),
    TP_Stage = unique(Stage[shortLetterCode == "TP"]),
    NT_Count = sum(shortLetterCode == "NT"),
    TP_Count = sum(shortLetterCode == "TP")
  ) %>%
  filter(!is.na(NT_Stage))

#### 결과 출력
print(stage_comparison)
cat("NT와 TP의 Stage가 일치하는 환자 수:", 
    sum(stage_comparison$NT_Stage == stage_comparison$TP_Stage, na.rm = TRUE), "\n")

duplicated_patients <- sum(duplicated(meta$patient))

cat("중복된 patient 개수:", duplicated_patients, "\n")

### duplicated
sum(duplicated(meta))

### Gender
gender <- table(meta$gender)
sum(is.na(meta$gender))
percent <- round(prop.table(gender) * 100, 2)  # prop.table로 비율 계산 후 퍼센트로 변환
labels_with_percent <- paste(names(gender), " (", percent, "%)", sep = "")
pie(gender, labels = labels_with_percent,  col = c("red", "blue"), main= "Gender")

### int data
for (i in int) {
  cat("=====", i, "=====\n")
  print(summary(meta[i]))
}

### 
row_df <- as.data.frame(row_data)
View(row_df)

# DEG analysis =================================================================
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

results2 <- topTable(fit2, coef = "Tumor_vs_Normal", adjust = "BH", number = Inf)
dim(results2)
dim(deg_results)

# Volcano Plot
results2$Significance <- ifelse(
  results2$adj.P.Val < 0.05 & abs(results2$logFC) > 1, "Significant", "Not significant"
  )

ggplot(results2, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("red", "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes (TP vs NT)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")



EnhancedVolcano(results2,
                lab = results2$ID,           # 유전자 이름
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
                widthConnectors = 0.5,             # 선 굵기
                max.overlaps = 5
)


head(rownames(deg_results))
head(rownames(exp_matrix_filtered))

all(rownames(exp_matrix_filtered) %in% rownames(deg_results))
all(rownames(deg_results) %in% rownames(exp_matrix_filtered))

head(rownames(v))


View(results2)

top_genes <- head(results2[order(results2$adj.P.Val), ], 30)$ID 
top_genes
top_exp <- exp_matrix_filtered[top_genes, ]

exp_long <- as.data.frame(t(top_exp))
exp_long$Sample <- rownames(exp_long)
exp_long <- reshape2::melt(exp_long, id.vars = "Sample")
colnames(exp_long) <- c("Sample", "Gene", "Expression")

exp_long$Group <- metadata_filtered$shortLetterCode[match(exp_long$Sample, rownames(metadata_filtered))]
View(exp_long)




Boxplot_gene <- function(df, g) {
  new_df <- subset(df, Gene == g)
  print(dim(new_df))
  p <- ggplot(new_df, aes(x = Gene, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +  # 이상치 표시 안함
    geom_jitter(position = position_jitter(0.2), alpha = 0.5) +  # 샘플 점 찍기
    theme_minimal() +
    scale_fill_manual(values = c("NT" = "blue", "TP" = "red")) +
    labs(title = "Expression of Top 30 DEGs (NT vs TP)",
         x = g,
         y = "Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  return(p)
}

exp_long$Expression <- log2(exp_long$Expression + 1)

Boxplot_gene(exp_long, top_genes[1])
Boxplot_gene(exp_long, top_genes[2])
Boxplot_gene(exp_long, top_genes[3])
Boxplot_gene(exp_long, top_genes[4])
Boxplot_gene(exp_long, top_genes[5])
Boxplot_gene(exp_long, top_genes[6])
Boxplot_gene(exp_long, top_genes[7])
Boxplot_gene(exp_long, top_genes[8])
Boxplot_gene(exp_long, top_genes[9])
Boxplot_gene(exp_long, top_genes[10])

# 
