---

---
___
# Differentially Expressed Genes analysis

특정 조건이나 실험처리에 따라 유전자 발현이 유의미하게 변화하는 유전자들을 식별하는 과정이다. 이 분석은 질병의 매커니즘을 이해하거나 새로운 치료 표적을 발견하는 데 중요한 역할을 한다.

### 원리
1. 정규화
	- 유전자 발현 데이터는 실험 간 또는 샘플 간의 비교가 가능하도록 정규화 과정을 거친다.
	- RNA-Seq 데이터는 각 샘플의 총 시퀀싱 뎁스나 유전자 길이 등을 고려하여 발현 수준을 조정
2. 통계적 분석
	- t-test, ANOVA 등의 방법을 사용해 각 유전자의 발현 수준이 통계적으로 유의미한 차이를 보이는지 평가합니다.
	- 일반적으로, 샘플 간의 발현 수준 차이가 통계적으로 유의미하다고 간주될 때 해당 유전자를 DEG로 분류합니다.
	- 다중 검정 문제를 고려하여, false discovery rate(FDR)를 조정합으로써, 잘못된 발견을 최소화합니다.
3. 생물학적 해석
	- 식별된 DEG들을 바탕으로, 이들이 참여하는 생물학적 과정이나 경로를 분석함으로써, 그 유전자들이 어떤 생물학적 기능을 수행하는지, 어떤 질병이나 생리적 조건과 관련이 있는지를 평가합니다.

### 주요 단계

GEOquery를 사용해 GEO(Gene Expression Omnibus)에서 유전자 발현 데이터셋을 다운로드하고 DEG 분석을 진행하려면, 데이터 수집부터 전처리, 분석, 시각화까지의 워크플로우를 설정해야 합니다. 아래는 R에서 GEOquery로 데이터를 가져와 DESeq2로 DEG 분석을 하고 그래프를 그리는 전체 과정을 단계별로 설명한 내용입니다.

---

### 1. **GEOquery로 데이터 다운로드**
GEO에서 데이터를 가져오려면 GEO accession 번호(예: GSE12345)를 알아야 합니다. 여기서는 예시로 "GSE157103"를 사용하겠습니다.

```R
# 필요한 패키지 설치 및 로드
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)

# GEO 데이터셋 다운로드 (예: GSE157103)
gse <- getGEO("GSE157103", GSEMatrix = TRUE, AnnotGPL = TRUE)

# 데이터셋 확인 (보통 리스트 형태로 반환됨)
gse <- gse[[1]]  # 첫 번째 데이터셋 선택 (GSE에 따라 다를 수 있음)
```

- **설명**: `getGEO`는 GEO 데이터셋을 가져오며, `GSEMatrix = TRUE`로 발현 매트릭스를 함께 다운로드합니다.

---

### 2. **데이터 전처리**
GEO 데이터는 발현 데이터와 메타데이터로 구성됩니다. 이를 DEG 분석에 맞게 준비합니다.

#### (1) 발현 데이터 추출
```R
# 발현 데이터 매트릭스 추출
exprs_data <- exprs(gse)  # 발현 값 (행: 유전자, 열: 샘플)

# 로그 변환 확인 (필요 시 변환)
if(max(exprs_data) > 50) {  # 로그 변환 안 된 경우
  exprs_data <- log2(exprs_data + 1)
}
```

#### (2) 메타데이터 추출
```R
# 샘플 메타데이터 추출
pheno_data <- pData(gse)

# 조건 열 확인 및 선택 (예: "condition" 열이 있다고 가정)
# 실제 데이터에 맞게 열 이름을 수정하세요
pheno_data$condition <- ifelse(grepl("normal", pheno_data$title, ignore.case = TRUE), "Normal", "GBM") 
colData <- data.frame(condition = pheno_data$condition, row.names = pheno_data$geo_accession) 
colData$condition <- factor(colData$condition, levels = c("Normal", "GBM")) # 예시 열 이름
# 데이터 정렬 
exprs_data <- exprs_data[, match(rownames(colData), colnames(exprs_data))]
```


#### (3) 데이터 정합성 확인
```R
# 샘플 이름이 발현 데이터와 메타데이터에서 일치하는지 확인
exprs_data <- exprs(gse) all(colnames(exprs_data) == rownames(colData)) # TRUE여야 함

all(colnames(exprs_data) == rownames(colData))  # TRUE여야 함
```

---

### 3. **DESeq2로 DEG 분석**
GEO 데이터가 RNA-Seq 카운트 데이터라면 DESeq2를 바로 사용할 수 있습니다. 하지만 마이크로어레이 데이터라면 `limma`를 사용하는 것이 더 적합할 수 있습니다. 여기서는 RNA-Seq 데이터라고 가정하고 진행합니다.

```R
# DESeq2 패키지 로드
library(DESeq2)

# DESeq2 데이터셋 생성
dds <- DESeqDataSetFromMatrix(countData = round(exprs_data),  # 정수형 카운트 필요
                              colData = colData,
                              design = ~ condition)

# 분석 실행
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", levels(colData$condition)[2], levels(colData$condition)[1]))

# DEG 필터링
deg <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
```

- **주의**: `exprs_data`가 이미 정규화된 값(예: 로그 변환된 마이크로어레이 데이터)이라면 DESeq2 대신 `limma`를 사용해야 합니다. RNA-Seq 원시 카운트 데이터일 경우만 DESeq2가 적합합니다.

---

### 4. **그래프 그리기**
DEG 분석 결과를 시각화합니다.

#### (1) Volcano Plot
```R
library(ggplot2)

# 결과 데이터프레임 변환
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "DEG", "Not DEG")

# Volcano Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("DEG" = "red", "Not DEG" = "grey")) +
  labs(title = "Volcano Plot of GSE157103", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
```

#### (2) Heatmap
```R
library(pheatmap)

# 정규화된 카운트 데이터 추출
norm_counts <- counts(dds, normalized = TRUE)
deg_counts <- norm_counts[rownames(deg), ]

# Heatmap 그리기
pheatmap(deg_counts, 
         scale = "row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Heatmap of DEGs in GSE157103",
         annotation_col = colData)
```

---

### 5. **실제 데이터에 맞게 조정**
- **GSE 확인**: 사용하려는 GSE 번호를 지정하세요(예: "GSE12345").
- **메타데이터 열 이름**: `pheno_data`에서 조건 열을 확인하고 `colData`에 맞게 수정하세요. 예를 들어, `pheno_data`를 출력해 열 이름을 확인할 수 있습니다:
  ```R
  head(pheno_data)
  ```
- **데이터 유형**: RNA-Seq이 아닌 마이크로어레이 데이터라면 아래와 같이 `limma`를 사용하세요:
  ```R
  library(limma)
  design <- model.matrix(~ condition, data = colData)
  fit <- lmFit(exprs_data, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, adjust = "fdr", number = Inf)
  deg <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
  ```

---

### 전체 코드 예시 (RNA-Seq 가정)
```R
# 패키지 로드
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# GEO 데이터 다운로드
gse <- getGEO("데이터셋이름", GSEMatrix = TRUE)[[1]]
exprs_data <- exprs(gse)
pheno_data <- pData(gse)
colData <- data.frame(condition = pheno_data$`treatment:ch1`)
colData$condition <- factor(colData$condition)

# DESeq2 분석
dds <- DESeqDataSetFromMatrix(countData = round(exprs_data), colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", levels(colData$condition)[2], levels(colData$condition)[1]))
deg <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Volcano Plot
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "DEG", "Not DEG")
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("DEG" = "red", "Not DEG" = "grey")) +
  labs(title = "Volcano Plot")

# Heatmap
norm_counts <- counts(dds, normalized = TRUE)
deg_counts <- norm_counts[rownames(deg), ]
pheatmap(deg_counts, scale = "row", main = "Heatmap of DEGs")
```

---

### 질문
- 사용하려는 GSE 번호가 있으면 알려주시면 그에 맞춰 조정할 수 있습니다.
- 데이터가 RNA-Seq인지, 마이크로어레이인지 확인되면 더 정확한 방법을 제안할게요! 추가로 궁금한 점 있으면 말씀해주세요.