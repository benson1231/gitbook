# Introduction of NGS

次世代定序（Next-Generation Sequencing, NGS）為現代基因體學研究提供了高通量、快速且成本效益高的 DNA/RNA 定序手段。以下概述其基本原理與主要定序技術類型。

---

## 一、NGS 基本原理

1. **文庫建置（Library Preparation）**：

   * 將 DNA/RNA 切割成小片段，接上接頭（adapters）。

2. **擴增（Amplification）**：

   * 使用 PCR 或橋式擴增進行片段數量擴大，以利訊號讀取。

3. **定序反應（Sequencing）**：

   * 利用酵素與螢光標記核苷酸讀取序列訊號。

4. **資料輸出（Base Calling）**：

   * 將影像訊號轉換為序列（FASTQ 格式）。

---

## 二、常見定序技術

| 技術平台                | 原理簡述                          | 特色             |
| ------------------- | ----------------------------- | -------------- |
| **Illumina**        | 合成測序（Sequencing by Synthesis） | 精準度高、讀長短、應用廣   |
| **Ion Torrent**     | 酸釋放偵測（半導體）                    | 成本低、反應快速       |
| **PacBio**          | 單分子即時測序（SMRT）                 | 可讀長達數萬 bp、誤差高  |
| **Oxford Nanopore** | 電流變化偵測 DNA 通過奈米孔時序列資訊         | 攜帶式設備、超長讀長、需校正 |

---

## 三、選擇依據與應用

* Illumina 適用於大多數短讀應用（exome、RNA-seq、ChIP-seq）
* PacBio 與 Nanopore 適用於全基因組組裝、結構變異分析
* Ion Torrent 適合目標區域定序、臨床檢測（如癌症 panel）

---

各種平台在精度、成本、讀長與資料處理上有所不同，研究者可根據研究目的（變異偵測、轉錄體分析、微生物分析等）選擇適合的定序技術，並搭配合適的生物資訊流程進行分析。
