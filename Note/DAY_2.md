# EDA 2
phenotype의 columns에 NAN을 어떻게 처리할 지 보기위해 NAN 값이 존재하는 columns 이며 데이터의 가치가 있는 columns들을 선별함
- colunms 조건
	1. value_count가 두가지 이상
	2. expression 데이터에 영향을 줄 수 있는 값
해당 조건 값들에 경우 대부분 NAN 값을 가지고 있다. 이에 column의 특성을 고려 해당 값의 NAN 값이 필연적으로 나올 수 있다는 것을 인지 해당 값을 대처할 수 있는 값은 대처 그중 가장 중요한 값은 `tumor grade:ch1`으로 암의 단계를 표현한 값이다. 하지만 해당 값은 데이터 타입이 문자열이기에 이를 수치로 환산하는 작업을 실행 또한 해당 값에 NAN 값을 해결할 수 있어 이를 진행

### **Tumor Grade (종양 등급) 개념 및 설명**

Tumor Grade(종양 등급)는 **암세포의 분화도(differentiation)와 증식 속도**를 기반으로 종양이 얼마나 공격적인지를 평가하는 지표이다. **일반적으로 숫자가 높을수록 악성도가 크고, 빠르게 증식하며 예후가 나쁘다**.

💡 **중요한 개념**

- **Low-grade (저등급, 분화도가 높음):** 정상 세포와 유사, 천천히 성장 (예: G1)
- **High-grade (고등급, 분화도가 낮음):** 비정상적이며 빠르게 성장 (예: G4)
- **Grade와 Stage의 차이:**
    - **Grade**: 암세포의 모양과 성장 속도에 따라 결정
    - **Stage**: 암이 신체에서 퍼진 정도를 나타냄 (예: TNM 단계)

---

## **1️⃣ 일반적인 Tumor Grade 분류 체계**

암의 종류에 따라 다르게 적용되지만, **세계보건기구(WHO)에서 정의한 4단계 시스템**이 가장 널리 쓰인다.

|Grade|설명|예후|
|---|---|---|
|**Grade I (G1)**|정상 조직과 유사한 구조, 천천히 성장, 낮은 악성도|예후 좋음|
|**Grade II (G2)**|약간 비정상적인 세포, 중간 정도 성장 속도|중간 수준|
|**Grade III (G3)**|정상 세포와 차이가 크고, 빠르게 증식|나쁜 예후|
|**Grade IV (G4)**|매우 비정상적, 공격적 성장, 침습성이 높음|매우 나쁜 예후|

---

## **2️⃣ 특정 종양에 따른 Grade 분류**

특정 암에서는 Grade 체계를 조금 다르게 적용하기도 한다.

### **🧠 뇌종양 (WHO Tumor Grading System)**

뇌종양(신경교종, glioma 등)에서는 WHO의 **Grade I ~ IV** 시스템을 사용함.

|WHO Grade|설명|예시|
|---|---|---|
|**Grade I (저등급, Low Grade)**|느린 성장, 비침습적|Pilocytic astrocytoma (소아에서 흔함)|
|**Grade II (중등급, Low to Intermediate Grade)**|불규칙한 세포 구조, 천천히 성장|Diffuse astrocytoma, OGII (oligodendroglioma II)|
|**Grade III (고등급, High Grade)**|빠른 증식, 악성도 증가|Anaplastic astrocytoma, OGIII|
|**Grade IV (고등급, Highly Malignant)**|빠른 성장, 높은 침습성|**Glioblastoma (GBM, 가장 공격적인 뇌종양)**|

💡 **예제에서 사용된 등급과 비교**

- **AII (Astrocytoma Grade II) → G2**
- **AIII (Astrocytoma Grade III) → G3**
- **AIII/IV → G3~G4 중간**
- **OGII (Oligodendroglioma Grade II) → G2**
- **OGIII (Oligodendroglioma Grade III) → G3**
- **IV (Glioblastoma 등) → G4**

📌 `AIII/IV` 같은 중간 단계는 예후와 성장을 고려하여 **G3.5 (또는 3)** 정도로 설정 가능.

---

## **3️⃣ 종양 분류 & 변환 가능성**

지금 제공한 등급을 **수치화**하면 아래와 같이 정리 가능하다.

|원본 등급|WHO Grade|변환 값 (수치)|
|---|---|---|
|Normal (NaN)|0|**0**|
|OGII|II (G2)|**1**|
|OGIII|III (G3)|**2**|
|AII|II (G2)|**1**|
|AIII|III (G3)|**2**|
|AIII/IV|III~~IV (G3~~G4)|**3** 또는 **2.5**|
|IV|IV (G4)|**4**|
|OA (Oligoastrocytoma)|II~~III (G2~~G3)|**2**|

👉 **즉, 저등급(low-grade)은 1, 고등급(high-grade)은 4로 변환 가능하며, 중간 등급은 보정 가능.**

---

## **4️⃣ 결론 및 최적의 변환 방식**

변환하는 기준을 **1~4 스케일**로 설정할 경우:

- **Normal = 0**
- **Grade II = 1**
- **Grade III = 2**
- **Grade III/IV = 3**
- **Grade IV = 4**

💡 **AIII/IV와 OA의 경우**

- `AIII/IV` = 2.5 또는 3
- `OA` = 2 또는 2.5

이렇게 설정하면 **연속적인 스케일로 분석하기 용이**하고, 머신러닝 모델이나 통계 분석에도 적절한 숫자 값으로 변환 가능! 🚀

NAN 값 중 cell line에 대한 값이 존재 이에
### **Cell Line에서 얻은 Tumor의 Grade (종양 등급) 지정 방법**

🔬 **Cell Line(세포주)에서 얻은 암세포는 일반적으로 "Grade"로 직접 분류되지 않는다.**  
그 이유는 **tumor grade는 조직학적 분류(환자의 실제 암 조직을 기반으로 하는 병리학적 평가)**를 따르기 때문이야. 하지만, **세포주의 유래(암 조직에서 유래된 단계)에 따라 암의 악성도 또는 분화도와 연결할 수 있다.**

---

## **1️⃣ 일반적인 Cell Line의 Tumor Grade 평가 방법**

세포주(Cell line)의 암 등급은 **보통 아래 세 가지 기준을 기반으로 유추**할 수 있어.

### **1) 원래 유래된 암 조직의 WHO Tumor Grade**

- 특정 **세포주가 환자의 어떤 암 조직에서 유래했는지**에 따라 대략적인 **암 등급(Grade I~IV)을 추정**할 수 있음.
- 예를 들어 **Glioblastoma Multiforme(GBM, 신경교종)**에서 유래된 세포주는 **WHO Grade IV**로 간주됨.
- 반면 **Oligodendroglioma Grade II에서 유래된 세포주는 WHO Grade II**로 평가 가능.

> **예시 (뇌종양 세포주)**
> 
> - **U87MG, LN229, T98G** → **Glioblastoma (WHO Grade IV)**
> - **SW1088** → **Astrocytoma (WHO Grade III)**
> - **HOG, MO3.13** → **Oligodendroglioma (WHO Grade II)**

➡ 즉, **세포주가 어떤 암에서 유래했는지가 중요**하며, 일반적으로 같은 종양 등급을 유지한다고 가정할 수 있음.

---

### **2) 암세포의 증식 속도 & 침습성**

**세포주의 증식 속도(Growth rate)와 침습성(Invasiveness)에 따라 대략적인 Grade를 추정 가능.**

- **고등급(High-grade, G3~G4) 종양**: 빠르게 증식하며, 공격적이고 침습성이 높은 세포주
    - 예: **HeLa (자궁경부암), MDA-MB-231 (삼중음성 유방암, TNBC), U87MG (뇌종양)**
- **저등급(Low-grade, G1~G2) 종양**: 분화도가 높고 천천히 자라는 세포주
    - 예: **MCF-7 (유방암, ER+), SW480 (대장암, 분화된 암종)**

**💡 예시:**

|Cell Line|유래|종양 Grade 추정|
|---|---|---|
|**MCF-7**|유방암 (ER+)|G2|
|**MDA-MB-231**|삼중음성 유방암 (TNBC)|G3~G4|
|**U87MG**|뇌종양 (Glioblastoma)|G4|
|**SW480**|대장암 (분화된 암종)|G2|
|**HCT116**|대장암 (미분화)|G3|

---

### **3) Cell Line의 유전자 돌연변이 및 표현형**

**암세포주의 유전자 변이 특징을 통해도 등급을 유추할 수 있어.**

- **P53, KRAS, EGFR 등의 돌연변이가 많고 불안정성이 높은 세포주** → 보통 **고등급(High-grade, G3~G4)**
- **분화 마커(예: CK18, CDH1 등)가 유지된 세포주** → 저등급(Low-grade, G1~G2)

예를 들어, **MDA-MB-231(유방암 세포주)**는 **P53 돌연변이가 존재하고, 증식이 빠르며, 침습성이 강하므로 G3~G4로 간주**됨.

---

## **2️⃣ 결론: Cell Line의 Tumor Grade 결정 방법**

1️⃣ **세포주가 유래한 암의 WHO Grade를 확인**

- Glioblastoma에서 유래 → **G4**
- Oligodendroglioma에서 유래 → **G2**

2️⃣ **증식 속도 & 침습성 분석**

- 빠른 증식 + 침습성 강함 → **G3~G4**
- 천천히 자라고 덜 침습적 → **G1~G2**

3️⃣ **유전자 돌연변이 정보 확인**

- P53/KRAS/EGFR 변이 多 → **G3~G4**
- 분화 마커 유지 → **G1~G2**

즉, **cell line 자체가 tumor grade를 직접 결정하는 것은 아니지만, 유래된 암 조직과 생물학적 특성을 바탕으로 grade를 추정할 수 있다.** 🧬💡

방법으로 처리

expression data와 phenotype data을 합치기 위해 expression의 row와 col을 전이함, 여기서 expression 데이터의 row 를 gene symbol로 바꾸었을 때 NAN로 되는 값이 이다는 것을 발견 해당 값들은 아직 gene symbol이 없는 값이나 실험에서 해당 실험이 잘 되었음을 보여 주는 데이터임을 확인 하지만 해당 데이터들을 없애기에 너무 많은 데이터로 사료되어 일단 원래 ID로 치환하기로 결정 

NAN 값은 ID로 유지하고 그렇기 않은 값은 gene symbol로 변경 후 전이 후 phenotype data 중 중요한 grade와 샘플이름을 합쳐 새로운 데이터프레임을 생성해 저장함