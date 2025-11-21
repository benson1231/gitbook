# TCGA

癌症基因體圖譜（The Cancer Genome Atlas, TCGA）是由美國國家癌症研究所（NCI）與國家人類基因組研究所（NHGRI）合作建立的大型癌症基因體資料庫。其目的是透過高通量技術（如 NGS）全面性分析多種癌症類型的分子變異，協助揭示癌症發展機制與臨床應用潛力。

---

## 1. TCGA 資料內容

TCGA 包含來自超過 30 種癌別、超過 11,000 位病人的資料，其涵蓋以下多種分子層級：

| 類型         | 說明                         |
| ---------- | -------------------------- |
| 基因體變異（DNA） | SNP、INDEL、CNV 等            |
| 表觀遺傳       | DNA 甲基化（450K / EPIC array） |
| 轉錄體（RNA）   | mRNA（RNA-seq）、miRNA        |
| 蛋白質層級      | RPPA 蛋白質定量平台               |
| 臨床資料       | 病患年齡、分期、預後、治療反應等           |

---

## 2. TCGA 資料取得與使用平台

### (1) GDC Portal

* 官方資料入口：[https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)
* 可瀏覽、篩選並下載所有 TCGA 資料（原始/標準化版本）

### (2) UCSC Xena

* [https://xenabrowser.net/](https://xenabrowser.net/)
* 支援多癌種整合瀏覽、即時分群與視覺化

### (3) R/Bioconductor 工具：`TCGAbiolinks`

* 可透過 R 程式碼自動化下載與分析 TCGA 資料

```r
library(TCGAbiolinks)

# 基因表現資料
query_exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
)
GDCdownload(query_exp)
exp_data <- GDCprepare(query_exp)

# 臨床資料
clinical_data <- GDCquery_clinic("TCGA-LUAD", "clinical")
```

---

## 3. 常見應用範例

| 分析類型           | 說明                            |
| -------------- | ----------------------------- |
| 差異表現分析（DEG）    | 比較癌症與正常組織、亞型之間的差異             |
| 存活分析（Survival） | 利用 Cox regression 分析基因與病人預後關聯 |
| 甲基化特徵分析        | 探索基因啟動子或特定區域的甲基化程度變化          |
| CNV / SNP 分析   | 分析基因擴增/刪除與點突變對癌症功能的影響         |
| 多組學整合分析        | 結合表現量、突變、臨床資料進行機器學習建模         |

---

## 4. 資料格式與範例

| 類別          | 格式                              | 範例檔名                        |
| ----------- | ------------------------------- | --------------------------- |
| Expression  | count / FPKM / TPM              | `HTSeq - Counts`            |
| Methylation | beta value matrix               | `methylation450k.tsv`       |
| Clinical    | JSON / TSV / XML                | `clinical_patient_brca.txt` |
| Mutation    | MAF（Mutation Annotation Format） | `*.maf.gz`                  |

---

## 5. 注意事項

* TCGA 多為美國族群資料，種族組成可能有偏差
* 某些癌別之正常對照組數量有限
* 臨床資料欄位名稱需經前處理標準化（如 vital\_status, days\_to\_death）
* 生物資訊分析中常搭配 GTEx 資料作為正常組背景

---

TCGA 是癌症基因體研究的基石，結合分子層級與臨床資訊，為癌症診斷、分類、預後與精準治療提供了寶貴資源。透過熟練操作其資料與工具，研究者能進行多維度的癌症系統生物學探索。
