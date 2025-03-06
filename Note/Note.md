# #1
딥러닝 모델을 어떻게 해서 만들 것인가?
- 현재 문제점
	- 데이터의 량이 생각보다 만다.
	- 이는 학습에


# matadata col
1. **barcode**: 샘플의 고유 식별자 (예: TCGA-XX-XXXX).
2. **patient**: 환자 ID.
3. **sample**: 샘플 ID (환자 내에서 여러 샘플이 있을 수 있음).
4. **shortLetterCode**: 샘플 유형을 나타내는 짧은 코드 (예: "TP" = Primary Tumor, "NT" = Normal Tissue).
5. **definition**: 샘플에 대한 설명 (예: "Primary solid Tumor").
6. **sample_submitter_id**: 제출자가 부여한 샘플 ID.
7. **sample_type_id**: 샘플 유형의 숫자 코드.
8. **tumor_descriptor**: 종양의 특성 (예: "Primary", "Metastatic").
9. **sample_id**: 샘플의 고유 ID.
10. **sample_type**: 샘플 유형 (예: "Tumor", "Normal").
11. **composition**: 샘플의 구성 (예: "Solid Tissue").
12. **days_to_collection**: 진단 후 샘플 채취까지의 일수.
13. **state**: 샘플 상태 (예: "released").
14. **initial_weight**: 샘플의 초기 무게.
15. **preservation_method**: 보존 방법 (예: "Frozen", "FFPE").
16. **intermediate_dimension**: 샘플 크기의 중간 치수.
17. **pathology_report_uuid**: 병리 보고서의 UUID.
18. **submitter_id**: 제출자 ID.
19. **shortest_dimension**: 샘플의 최소 치수.
20. **oct_embedded**: OCT(Optimal Cutting Temperature) 매체 사용 여부.
21. **specimen_type**: 검체 유형.
22. **longest_dimension**: 샘플의 최대 치수.
23. **is_ffpe**: FFPE(Formalin-Fixed Paraffin-Embedded) 여부.
24. **tissue_type**: 조직 유형 (예: "Tumor", "Normal").
25. **synchronous_malignancy**: 동시 악성 종양 여부.
26. **ajcc_pathologic_stage**: AJCC 병리학적 병기 (예: "Stage IA").
27. **days_to_diagnosis**: 출생 후 진단까지의 일수.
28. **laterality**: 병변의 위치 (예: "Left", "Right").
29. **treatments**: 치료 정보.
30. **tissue_or_organ_of_origin**: 종양의 기원 조직/장기.
31. **age_at_diagnosis**: 진단 시 나이.
32. **primary_diagnosis**: 주 진단명 (예: "Lung Adenocarcinoma").
33. **prior_malignancy**: 이전 악성 종양 여부.
34. **year_of_diagnosis**: 진단 연도.
35. **prior_treatment**: 이전 치료 여부.
36. **diagnosis_is_primary_disease**: 진단이 주요 질병인지 여부.
37. **ajcc_staging_system_edition**: AJCC 병기 시스템 버전.
38. **ajcc_pathologic_t**: AJCC T 병기 (종양 크기/침범 정도).
39. **morphology**: 종양의 형태학적 분류 (ICD-O 코드).
40. **ajcc_pathologic_n**: AJCC N 병기 (림프절 전이).
41. **ajcc_pathologic_m**: AJCC M 병기 (원격 전이).
42. **residual_disease**: 잔류 질병 상태.
43. **classification_of_tumor**: 종양 분류.
44. **diagnosis_id**: 진단 ID.
45. **icd_10_code**: ICD-10 질병 코드.
46. **site_of_resection_or_biopsy**: 절제/생검 부위.
47. **tumor_of_origin**: 종양 기원.
48. **sites_of_involvement**: 침범 부위.
49. **tobacco_smoking_quit_year**: 흡연 중단 연도.
50. **tobacco_smoking_status**: 흡연 상태 (예: "Current", "Never").
51. **exposure_id**: 노출 정보 ID.
52. **exposure_type**: 노출 유형.
53. **pack_years_smoked**: 흡연량 (팩-년 단위).
54. **tobacco_smoking_onset_year**: 흡연 시작 연도.
55. **race**: 인종.
56. **gender**: 성별.
57. **ethnicity**: 민족.
58. **vital_status**: 생존 상태 (예: "Alive", "Dead").
59. **age_at_index**: 기준 시점에서의 나이.
60. **days_to_birth**: 출생일로부터의 일수 (음수로 표현됨).
61. **demographic_id**: 인구통계 ID.
62. **age_is_obfuscated**: 나이 정보가 모호화되었는지 여부.
63. **country_of_residence_at_enrollment**: 등록 시 거주 국가.
64. **days_to_death**: 사망까지의 일수.
65. **bcr_patient_barcode**: BCR(Biomedical Research) 환자 바코드.
66. **primary_site**: 주요 부위.
67. **project_id**: 프로젝트 ID (예: "TCGA-LUAD").
68. **disease_type**: 질병 유형.
69. **name**: 이름 (보통 프로젝트/데이터셋 이름).
70. **releasable**: 데이터 공개 가능 여부.
71. **released**: 데이터 공개 여부.
72. **figo_stage**: FIGO 병기 (주로 여성 암에 사용).
73. **figo_staging_edition_year**: FIGO 병기 버전 연도.
74. **paper_patient**: 논문에서 사용된 환자 ID.
75. **paper_Sex**: 논문에서의 성별.
76. **paper_Age.at.diagnosis**: 논문에서의 진단 시 나이.
77. **paper_T.stage**: 논문에서의 T 병기.
78. **paper_N.stage**: 논문에서의 N 병기.
79. **paper_Tumor.stage**: 논문에서의 종양 병기.
80. **paper_Smoking.Status**: 논문에서의 흡연 상태.
81. **paper_Survival**: 논문에서의 생존 정보.
82. **paper_Transversion.High.Low**: 돌연변이 패턴 (High/Low).
83. **paper_Nonsilent.Mutations**: 비침묵 돌연변이 수.
84. **paper_Nonsilent.Mutations.per.Mb**: Mb당 비침묵 돌연변이 수.
85. **paper_Oncogene.Negative.or.Positive.Groups**: 온코진 그룹.
86. **paper_Fusions**: 융합 유전자 정보.
87. **paper_expression_subtype**: 발현 아형.
88. **paper_chromosome.affected.by.chromothripsis**: 염색체 파괴 현상.
89. **paper_iCluster.Group**: iCluster 그룹.
90. **paper_CIMP.methylation.signature**: CIMP 메틸화 서명.
91. **paper_MTOR.mechanism.of.mTOR.pathway.activation**: mTOR 경로 활성화 메커니즘.
92. **paper_Ploidy.ABSOLUTE.calls**: 배수성(ABSOLUTE 호출).
93. **paper_Purity.ABSOLUTE.calls**: 순도(ABSOLUTE 호출).
94. **Stage**: 병기 (사용자 정의 열일 가능성 있음).


# 결측치 처리
### 2. Stage의 NA 처리 가능성

Stage가 NA인데 TNM 데이터가 존재하므로, 이를 활용해 병기를 추정할 수 있습니다. AJCC 병기 분류 기준에 따라 T, N, M 값을 조합해 Stage를 예측할 수 있습니다.

#### AJCC 폐암 병기 기준 (간략화)

- **Stage I**: T1 (T1a, T1b, T1c) + N0 + M0.
- **Stage II**: T2 (T2a, T2b) + N0 + M0 또는 T1/T2 + N1 + M0.
- **Stage III**: T3/T4 + N0/N1 + M0 또는 T1~T4 + N2/N3 + M0.
- **Stage IV**: M1 (원격 전이).

#### 데이터 기반 추정

- **T1a (2개)**: N0, M0라면 Stage IA.
- **T1c (1개)**: N0, M0라면 Stage IA.
- **T2 (1개), T2a (1개), T2b (2개)**: N0, M0라면 Stage IB 또는 IIA.
- **T3 (1개)**: N0, M0라면 Stage IIB.
- **T4 (1개)**: N0, M0라면 Stage IIIA 가능.
- **M1 (1개)**: Stage IV 확정.

