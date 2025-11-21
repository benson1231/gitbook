# Bioinformatics Tools for NGS

次世代定序（Next-Generation Sequencing, NGS）產生大量高通量序列資料，需依賴一系列生物資訊工具進行前處理、比對、變異分析與功能註解。以下介紹常見工具及其應用階段。

---

## 一、常見分析流程概覽

1. **資料品質檢查（QC）**

   * FastQC：檢視 read 品質、GC 分布、adapter contamination 等
   * MultiQC：整合多個樣本 QC 報告

2. **序列修剪與過濾**

   * Trimmomatic / Cutadapt / fastp：去除低品質序列與接頭污染

3. **比對至參考基因組**

   * BWA / Bowtie2：短讀序列比對（DNA-seq）
   * STAR / HISAT2：轉錄體資料比對（RNA-seq）

4. **SAM/BAM 處理**

   * SAMtools / Picard：格式轉換、排序、去除 duplicates

5. **變異偵測（Variant Calling）**

   * GATK（HaplotypeCaller）：常用於 germline variant 分析
   * FreeBayes / DeepVariant：替代工具
   * VarScan2 / Mutect2：腫瘤樣本 somatic variant 分析

6. **註解與解釋**

   * ANNOVAR / VEP：將 SNP/INDEL 加入基因與疾病資訊
   * SnpEff：功能影響預測

---

## 二、轉錄體分析（RNA-seq）

* FeatureCounts / HTSeq：定量基因表現
* DESeq2 / edgeR：差異表現分析
* clusterProfiler / GSEA：路徑富集分析

---

## 三、表觀遺傳 / ChIP-seq / ATAC-seq

* MACS2：peaks 偵測
* DeepTools：coverage、可視化
* ChIPseeker：結合註解分析

---

## 四、多樣性與菌叢分析（16S/Metagenomics）

* QIIME2 / mothur：微生物族群結構分析
* Kraken2 / MetaPhlAn：物種分類與豐度估計
* HUMAnN3：功能路徑分析

---

## 五、常見輔助工具與平台

| 工具／平台              | 用途說明                     |
| ------------------ | ------------------------ |
| Galaxy             | 圖形化介面分析平台，適合非程式背景用戶      |
| Bioconda           | 安裝 NGS 工具的 Conda 套件管理資源庫 |
| Nextflow/Snakemake | 建立重複性與模組化的分析流程           |
| IGV                | 基因組視覺化瀏覽器                |

---

這些工具組成一套完整的 NGS 生物資訊分析流程，可依據應用目的（癌症研究、轉錄分析、病原檢測等）做適當組合與最佳化，是當代基因體研究的核心能力。
