# DEG
# library insall 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("DEGseq")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")

install.packages("pheatmap")

library(GEOquery)
library(tibble)
library(DEGseq)
library(ggplot2)
library(dplyr)
library(splines)
library(limma)
library(EnhancedVolcano)
# 필요한 패키지 로드
set.seed(123)  # 재현성 보장


# 1. GEO 데이터 다운로드 및 선택
gset <- getGEO("GSE15824", GSEMatrix = TRUE, AnnotGPL = TRUE)
idx <- grep("GPL570", names(gset), ignore.case = TRUE)  # GPL570 플랫폼 선택
if (length(idx) == 0) idx <- 1  # 플랫폼 없을 경우 첫 번째 데이터셋 사용
gset <- gset[[idx]]

# 열 이름 중복 방지 (fvarLabels 정제)
fvarLabels(gset) <- make.names(fvarLabels(gset))

# 2. 샘플 그룹 정의 및 필터링
gsms <- "000000000000XXXXXXXXXXXXXXXXXXXXX11XXXXXXXXXX"  # 질문에서 제공된 샘플 마스크
sml <- strsplit(gsms, "")[[1]]
sel <- which(sml != "X")  # "X" 제외한 샘플 선택
gset <- gset[, sel]  # 선택된 샘플만 남김
group_labels <- ifelse(sml[sel] == "0", "GBM", "1")  # "0" = GBM, "1" = normal
group_labels[group_labels == "1"] <- "normal"  # "1"을 "normal"로 변환

# 3. 발현 데이터 전처리
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)  # 로그 변환 필요 여부 판단
if (LogC) {
  ex[ex <= 0] <- NaN  # 0 이하 값을 NaN으로 변환
  exprs(gset) <- log2(ex)  # 로그2 변환 적용
}

# 4. 그룹 설정 및 디자인 매트릭스
gs <- factor(group_labels, levels = c("normal", "GBM"))  # "normal"을 기준으로 설정
pData(gset)$group <- gs  # gset 객체에 그룹 정보 추가
design <- model.matrix(~ group + 0, data = pData(gset))  # Intercept 없는 모델
colnames(design) <- levels(gs)  # 열 이름 설정: "normal", "GBM"

# 5. limma를 사용한 DEG 분석
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts("GBM-normal", levels = design)  # GBM vs normal 대조
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend = TRUE)  # 트렌드 반영으로 robust 결과 도출

# 6. DEG 결과 추출
deg_results <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
deg_results$Gene.symbol <- fData(gset)$`Gene symbol`[match(rownames(deg_results), 
                                                           rownames(fData(gset)))]  # 유전자 심볼 추가
# fData 열 이름 확인 및 유전자 심볼 추가
colnames(fData(gset))  # 열 이름 출력하여 확인
# "Gene symbol"이 없으면 적절한 열 이름으로 대체 (예: "SYMBOL" 또는 "Gene.Symbol")
gene_col <- "Gene symbol"  # 기본값
if (!("Gene symbol" %in% colnames(fData(gset)))) {
  gene_col <- grep("symbol", colnames(fData(gset)), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(gene_col)) gene_col <- rownames(fData(gset))  # 심볼 없으면 probe ID 사용
}
deg_results$Gene.symbol <- fData(gset)[[gene_col]][match(rownames(deg_results), rownames(fData(gset)))]

# 7. EnhancedVolcano 플롯 생성
EnhancedVolcano(deg_results,
                lab = deg_results$Gene.symbol,         
                x = 'logFC',                          
                y = 'adj.P.Val',                      
                pCutoff = 0.05,                       
                FCcutoff = 1,                         
                title = 'Differential Expression: GBM vs Normal (GSE15824)',
                subtitle = 'Microarray Data Analysis',
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 Adjusted P-value',
                pointSize = 2,
                labSize = 3,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendPosition = 'top',
                legendLabels = c('NS', 'Log2FC > 1', 'FDR < 0.05', 'FDR < 0.05 & Log2FC > 1'),
                selectLab = c("FN1", "CD44", "MYC", "CDK1", "SERPINE1", 
                              "COL3A1", "COL1A2", "LOX", "POSTN", "EZH2"))

# 8. 결과 저장
write.csv(deg_results, "GSE15824_DEG_results.csv", row.names = FALSE)
ggsave("GSE15824_volcano_plot.png", width = 10, height = 8, dpi = 300)