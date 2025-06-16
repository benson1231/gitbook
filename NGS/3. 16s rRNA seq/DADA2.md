# DADA2

DADA2（Divisive Amplicon Denoising Algorithm）是一套用於 16S rRNA 與其他擴增子定序資料的錯誤修正與序列推斷工具，核心理念是：

> 「不依賴參考資料庫，透過建模錯誤機率，直接從定序資料中推斷真實存在的變異序列（ASV）」

DADA2 最大特色是能辨識到單一鹼基變異（single nucleotide resolution），產出精確的 Amplicon Sequence Variants（ASVs），而非粗略的 OTU 聚類。

---

### 🔬 DADA2 核心演算法流程

1. **品質過濾與截尾（filterAndTrim）**

   * 去除低品質 reads、裁剪 primers/adapters。

2. **學習錯誤模型（learnErrors）**

   * 根據輸入資料計算不同鹼基在不同位置上的誤讀機率，建立錯誤矩陣。

3. **去噪（dada）**

   * 利用錯誤模型推斷每條讀段是否為真實序列或為錯誤產物。

4. **配對合併（mergePairs）**

   * 將 forward 與 reverse reads 對齊合併，並去除不一致區段。

5. **去除嵌合體（removeBimeraDenovo）**

   * 偵測並移除因 PCR 錯配造成的嵌合序列。

6. **建立 ASV 表（makeSequenceTable）**

   * 得到樣本 x ASV 的 abundance 矩陣。

7. **分類註解（assignTaxonomy）**

   * 使用參考資料庫（如 SILVA）比對每條 ASV，取得分類資訊。

---

### ✅ DADA2 相較傳統 OTU 聚類的優勢

| 面向       | OTU 方法     | DADA2 / ASV 方法 |
| -------- | ---------- | -------------- |
| 分群原則     | 序列相似度 ≥97% | 單一鹼基解析度        |
| 是否依賴參考序列 | 是/否皆可      | 不需要（只在註解階段才需）  |
| 可重現性     | 不穩定（依聚類參數） | 高（序列為唯一單位）     |
| 結果解釋力    | 低（混合菌群）    | 高（可精確追蹤變異）     |

---

## 🧬 DADA2 實際應用於 16S rRNA 分析

參見[DADA技術文件](https://benjjneb.github.io/dada2/index.html)

### 最終可得到：

* `seqtab.nochim`：ASV abundance 表格（樣本 × 序列）
* `taxa`：每個 ASV 對應的分類名稱（Kingdom → Genus）

這些結果可輸入至 `phyloseq` 進行後續多樣性分析、可視化、統計比較等。

---

接下來請參考本文件其他章節來進行 ASV 表與分類結果的多樣性分析與視覺化。
