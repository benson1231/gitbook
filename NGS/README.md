# 🧬 Next-Generation Sequencing

次世代定序（Next-Generation Sequencing, NGS）是一種高通量 DNA/RNA 序列定序技術，能在短時間內產出大量資料，應用於基因體學、轉錄體學、腫瘤分析、微生物定序等多個領域。

---

## 1. NGS 與傳統定序的差異

| 項目    | 傳統 Sanger 定序 | NGS（次世代定序）     |
| ----- | ------------ | -------------- |
| 產出量   | 單一序列         | 每次數百萬條序列       |
| 成本與時間 | 高            | 相對低廉且快速        |
| 應用範圍  | 基因驗證、小規模分析   | 全基因組、轉錄體、大規模分析 |

---

## 2. NGS 工作流程

### (1) 样本準備

* 提取 DNA 或 RNA
* 建立定序資料庫（library）：加入接頭（adapters）、酶切或 PCR 擴增

### (2) 定序（Sequencing）

* 使用平台如 Illumina、Ion Torrent、Nanopore
* 根據技術不同產生短讀（short reads）或長讀（long reads）

### (3) 資料處理（Bioinformatics）

* **質控（Quality Control）**：如 FastQC
* **比對（Alignment）**：如 BWA、STAR、HISAT2
* **變異檢測（Variant Calling）**：如 GATK、FreeBayes
* **表現量分析（Expression Quantification）**：如 featureCounts、Salmon
* **可視化與統計分析**：如 IGV、R/Seurat

---

## 3. 常見 NGS 應用

| 應用類型                          | 說明                      |
| ----------------------------- | ----------------------- |
| Whole Genome Sequencing (WGS) | 完整基因體定序                 |
| Whole Exome Sequencing (WES)  | 只定序編碼區（外顯子）             |
| RNA-seq                       | 分析基因表現量與轉錄變異            |
| Small RNA-seq                 | 分析 miRNA、siRNA 等小分子 RNA |
| Single-cell RNA-seq           | 單細胞層級分析轉錄體              |
| ChIP-seq                      | 分析蛋白質與 DNA 的結合區域（如轉錄因子） |

---

## 4. NGS 平台比較

| 平台              | 讀長類型 | 優點             | 缺點       |
| --------------- | ---- | -------------- | -------- |
| Illumina        | 短讀   | 高準確率、高通量       | 難以分析重複區域 |
| Oxford Nanopore | 長讀   | 可即時分析、便攜裝置     | 錯誤率較高    |
| PacBio          | 長讀   | 讀長穩定、適合全長轉錄本分析 | 成本較高     |

---

## 5. 資料儲存格式

| 格式      | 說明                      |
| ------- | ----------------------- |
| FASTQ   | 原始讀取資料，包含序列與品質分數        |
| BAM/SAM | 對齊結果，SAM 為文字，BAM 為二進位格式 |
| VCF     | 變異資料格式（如 SNP、INDEL）     |
| GTF/GFF | 基因註解格式                  |

---

NGS 技術已成為現代生物醫學研究與精準醫療的重要支柱。掌握其流程與應用對於從事基因體與轉錄體研究至關重要。
