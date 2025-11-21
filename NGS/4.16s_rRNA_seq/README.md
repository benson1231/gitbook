# Full-Length 16S rRNA Sequencing Analysis Pipeline

本文件說明從原始資料到進階統計分析的完整 16S rRNA 全長定序資料分析流程，並提供各步驟的簡要說明。

---

## 1. 原始資料處理（Raw Data Processing）

### ASV 產生（ASVs Generation）

* **品質過濾（Quality Filtering）**：去除低品質讀段、修剪接頭與引子，以提升後續分析的準確度。
* **去噪（Denoising，使用 DADA2）**：辨識並修正定序錯誤，產生唯一的 Amplicon Sequence Variants（ASVs），精確度高於傳統 OTUs。

### 註解（Annotation）

* 將 ASV 定序比對至微生物資料庫，進行分類學註解。

  * **NCBI**：通用參考庫。
  * **GreenGenes**：較舊但廣泛使用的分類庫。
  * **SILVA**：更新頻繁、分類精細，廣泛應用。
  * **HOMD**：人類口腔微生物專用資料庫。
  * **UNITE**：真菌 ITS 專用資料庫。

**輸出**：ASV 表格，包含每個樣本中出現的 ASVs 及其分類資訊。

---

## 2. 下游分析模組（Downstream Analysis Modules）

### Alpha 多樣性（Alpha Diversity，樣本內多樣性）

評估單一樣本中菌相的豐富度與均勻度。

* **多樣性指數**：如 Shannon、Simpson、Chao1 等指標。
* **稀疏曲線（Rarefaction Curve）**：評估測序深度與物種數之關係。
* **豐度排序圖（Rank Abundance Curve）**：視覺化菌相分布。
* **物種累積曲線（Specaccum Curve）**：反映發現新物種的速率。
* **Shannon-Wiener 曲線**：綜合反映豐富度與均勻性。

### Beta 多樣性（Beta Diversity，樣本間差異）

比較樣本之間的菌相差異與群落結構。

* **PCA / 3D PCA**：線性降維分析。
* **PCoA / 3D PCoA**：基於距離矩陣（如 UniFrac）。
* **NMDS**：非參數排序方法。
* **t-SNE**：非線性降維，適用於高維資料。
* **PLS-DA**：監督式分群分析。
* **UPGMA Tree**：分層聚類法。
* **UniFrac 距離盒狀圖**：顯示演化關係差異。

### 系統發生分析（Phylogenetic Analysis）

* **演化樹構建（Phylogenetic Tree）**：利用比對後序列建立菌株間的系統發生樹。

### 群落組成分析（Community Composition）

視覺化菌群豐度與結構。

* **Venn / Upset Plot**：比較組間共享與專一菌種。
* **Heatmap**：菌種相對豐度矩陣圖。
* **分類組成（Taxonomy Profile）**：以柱狀圖呈現不同分類層級（門、屬等）。
* **KRONA**：互動式圓形分類圖。
* **TreeMap / HeatTree**：以區域大小或樹狀結構呈現豐度。
* **Ternary Plot**：三組樣本之間菌群差異視覺化。

### 關聯分析（Correlation Analysis）

分析菌種之間的共現或相互作用。

* **Spearman 關聯分析**：非參數相關檢定。
* **網絡分析（Network Analysis）**：建立菌種間關聯圖譜。
* **相關熱圖（Correlation Heatmap）**：顯示菌種間關係。

### 功能預測（Functional Prediction）

預測微生物群落可能的功能。

* **PICRUSt2**：基於 16S 資料推測代謝途徑。
* **Tax4Fun**：整合 SILVA 與 KEGG 進行功能注釋。
* **FAPROTAX**：針對生地化功能進行預測。

### 統計分析（Statistical Analysis）

進行統計檢定與生物標記探索。

* **Kruskal-Wallis**：非參數多組檢定。
* **STAMP**：簡易圖形化統計分析工具。
* **ANCOM**：針對組成性資料的差異分析。
* **ALDEx2**：貝式差異豐度分析法。
* **LEfSe**：LDA 為基礎的生物標誌物篩選方法。
* **MetagenomeSeq**：處理稀疏計數資料。
* **Bubble Chart / Dominant Taxa**：視覺化主要菌種變化。
* **ANOSIM / MRPP / Adonis**：多變量群組差異分析。

### 環境因子分析（Environment Factor Analysis）

探討環境或臨床變數與菌相的關聯。

* **RDA / CCA / dbRDA**：約束排序方法。
* **Mantel Test**：距離矩陣相關性分析。
* **BioENV**：篩選解釋變異的重要變項。

### 進階分析（Advanced Analysis）

* **MicroPITA**：找出具區辨力的核心菌種。
* **Source Tracker**：推測微生物來源（如糞便、皮膚等）。

---

## 工具建議（Tools & Recommendations）

* **QIIME2**：完整工作流程工具，適合全流程分析。
* **R 語言套件**：如 `phyloseq`、`vegan`、`microbiome` 等，適合彈性視覺化與統計分析。
* **Galaxy 平台**：免程式操作，適合初學者。
* **模擬樣本**：ZymoBIOMICS 等提供可重複測試用菌群，適合方法驗證。

---

## 附註（Notes）

* 實驗設計應包含陰性與陽性對照。
* 分析時務必結合 metadata 資訊。
* 多次重複樣本有助於提高統計可靠性。

---

若需練習數據，可參考 NCBI SRA、Qiita 資料庫，或向我們請求針對特定微生物群（如腸道、口腔、環境）模擬產生的 FASTQ 數據集。
