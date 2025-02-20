# Glioblastoma Multiforme, GBM
-> GBM, 교모세포종은 성인의 노에서 가장 흔하고 공격적인 원발성악성 신경교종이다.

---
### 특징
- 악성 등급 : WHO 4등급 (가장 높은 등급의 신경교종)
- 발생 부위 : 대뇌 반구에서 주로 발생하지만, 뇌의 다른 부분으로 확산될 수 있음
- 성장 속도 : 빠르고 침습적인 성장 양상을 보이며, 주변 정장 뇌조직과 명확한 경계를 가지지 않음
- 생존율 : 진단 후 평균 생존 기간이 12 ~ 15 개월로 매우 낮음

---
### 원인 및 발생 기전
- 유전적 요인 : TP53, PTEN, EGFR, IDH1/2
- 분자적 특징 :
	- IDH 변이 여부 : IDH 변이가 없느 GBM은 더욱 공격적인 경향
	- EGFR 증폭 : 종양 성장 촉진
	- MGMT 프로모터 메틸화 : 항암제에 대한 반응ㅇ성 증가
- 환경적 요인 : 방사선 노출, 유전적 소인 등이 영향을 줄 수 있음

---

# Materials and Methods
### Install package and load packages
1️⃣`library(DEGseq)`
- 차등 유전자 발현 분석(Differential Gene Expression, DEG) 수행
- RNA-Seq 데이터에서 유의한 발현 차이를 보이는 유전자를 찾는 데 사용

2️⃣ `library(GEOquery)`
- NCBI의 **GEO (Gene Expression Omnibus)** 데이터베이스에서 마이크로어레이 및 RNA-Seq 데이터를 다운로드하는 데 사용
- `getGEO()` 함수를 이용해 원하는 데이터셋을 가져올 수 있음

3️⃣ `library(ggplot2)`
- 데이터 시각화 라이브러리
- DEG 분석 후 Volcano Plot, MA Plot 등을 그릴 때 유용

4️⃣ `library(dplyr)`
- 데이터 전처리 및 조작을 위한 필수 패키지
- `filter()`, `mutate()`, `group_by()`, `summarize()` 등의 함수 활용 가능

5️⃣ `library(splines)`
- **스플라인(Spline) 보간법 및 회귀 모델링**을 위한 패키지
- 비선형 데이터 분석 및 시각화에 사용됨

``` R
library(GEOquery)
library(DEGseq)
library(ggplot2)
library(dplyr)
library(splines)
```

###  Data load
``` R
gse <- getGEO("GSE15824") # 데이터를 가져옴
pheno <- attr(attr(ges[[1]], "phenoData"), "data") # 표현형 데이터를 추출
exprs <- attr(ges[[1]], "assayData")[["exprs"]] # 유전자 발현 데이터를 추출
```
#####  GEO 데이터
GEO (Gene Expression Omnibus)에서 **GSE**와 **GPL**은 각각 **실험 데이터**와 **플랫폼(기술)**을 의미합니다.
##### 🔹 GSE (GEO Series)

- **GSE (GEO Series Expression)**는 **하나의 연구(실험)에서 생성된 데이터셋**을 의미합니다.
- 보통 하나의 **GSE**에는 여러 샘플(환자군, 대조군 등)이 포함됩니다.
- 연구자가 실험을 수행한 후 GEO에 제출하는 **데이터 단위**입니다.

📌 **예제**:

- `GSE183947` → 특정 연구에서 수행된 **마이크로어레이/RNA-Seq 실험 데이터셋**
- `GSEXXXXX` 형태로 명명되며, 같은 실험에서 얻어진 여러 샘플이 포함됨.

📌 **GSE에 포함된 정보**:

1. **실험 조건** (예: 대조군 vs 환자군)
2. **샘플 ID (GSM) 목록**
3. **사용한 기술 (GPL)**
4. **유전자 발현 값 (expression matrix)**
5. **논문 정보 (관련 연구 논문 링크 포함)**

---

##### 🔹 GPL (GEO Platform)

- **GPL (GEO Platform)**은 **사용된 마이크로어레이 칩 또는 시퀀싱 기술**을 의미합니다.
- 연구자가 어떤 기술로 데이터를 생성했는지 나타내는 **플랫폼 ID**입니다.

📌 **예제**:

- `GPL570` → **Affymetrix Human Genome U133 Plus 2.0 Array**
- `GPL10558` → **Illumina HumanHT-12 V4.0 expression beadchip**
- `GPLXXXX` 형태로 명명되며, 같은 기술을 사용한 여러 연구에서 공유됨.

📌 **GPL에 포함된 정보**:

6. **플랫폼 제조사 (Affymetrix, Illumina 등)**
7. **프로브 ID 및 유전자 매핑 정보**
8. **기술 유형 (Microarray, RNA-Seq 등)**
9. **참고 논문 및 기술 문서**

---

##### **📌 GSE와 GPL 관계**

- **GSE는 연구(실험) 단위**이고, **GPL은 해당 연구에서 사용한 기술(플랫폼)**입니다.
- 하나의 **GPL(플랫폼)**은 여러 **GSE(연구 데이터셋)**에서 사용될 수 있습니다.

📌 **예제 (관계 이해하기)**

10. `GSE183947` → 유방암 관련 연구 데이터셋
    
    - **사용된 플랫폼**: `GPL570` (Affymetrix 마이크로어레이 칩)
11. `GSEXXXXX` → 다른 연구 데이터셋
    
    - **사용된 플랫폼**: `GPL570` (같은 플랫폼 사용)
    - 즉, **다른 연구지만 같은 기술을 사용했을 수 있음**.

---

##### **📌 GEO에서 데이터를 검색하는 방법**

12. **GSE 검색 (특정 연구 찾기)**
    
    ```r
    library(GEOquery)
    gse <- getGEO("GSE183947", GSEMatrix = TRUE)
    show(gse)
    ```
    
13. **GPL 검색 (플랫폼 정보 확인)**
    
    ```r
    gpl <- getGEO("GPL570")
    show(gpl)
    ```
    

---

##### **🔬 연구 적용**

- **"GSE 찾기"** → 특정 암 연구 데이터를 찾고 싶을 때 사용
- **"GPL 확인"** → 데이터의 유전자 발현 값이 어떤 기술에서 나왔는지 확인할 때 사용
- **"GSE-GPL 매핑"** → 같은 플랫폼을 사용한 다른 연구와 비교 분석할 때 유용

🚀 **즉, GSE는 연구 데이터**, **GPL은 데이터 생성 기술(플랫폼)**입니다!

- GSE15824
	- 30 개
		- GBM : 12개
		- normal : 5개
