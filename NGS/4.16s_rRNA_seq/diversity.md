# Diversity

微生物群落的多樣性是環境與健康研究中極具價值的指標。透過 16S rRNA 定序技術，我們可以深入了解樣本內（Alpha 多樣性）與樣本間（Beta 多樣性）的菌群組成差異。不同的多樣性指標與統計方法，能協助我們評估微生物群落的豐富程度、均勻性，以及群落結構的穩定性與變化。

---

## 🧬 Alpha 多樣性分析（Alpha Diversity）

Alpha 多樣性是指「**單一樣本內的菌群多樣性**」，通常用於衡量樣本中的微生物豐富度與均勻度。

### 常見指標：

* **Observed OTUs / ASVs**：實際觀測到的物種數量。
* **Chao1**：估計物種豐富度，強調低豐度物種。
* **Shannon index**：考慮物種數量與分布均勻程度。
* **Simpson index**：衡量優勢物種的主導程度。

### 可視化方式：

* **稀疏曲線（Rarefaction curve）**：觀察物種數是否已飽和。
* **箱型圖（Boxplot）**：比較不同組別樣本之間的多樣性差異。

### 目的：

* 瞭解樣本本身菌相是否豐富、多樣、或由單一物種主導。
* 可用於比較 treatment/control、疾病/健康組間之差異。

---

## 🧭 Beta 多樣性分析（Beta Diversity）

Beta 多樣性是指「**不同樣本之間的菌群組成差異**」，用於評估整體菌群結構的變異。

### 距離計算方法：

* **Bray-Curtis 距離**：考慮物種豐度的相異度。
* **Jaccard 距離**：僅考慮物種有無。
* **UniFrac（weighted/unweighted）**：結合系統發生關係的距離（基於演化樹）。

### 常見降維方法與圖示：

* **PCoA（主座標分析）**
* **NMDS（非度量多維尺度分析）**
* **PCA（主成分分析）**
* **t-SNE / UMAP（高維視覺化）**

### 群組比較：

* **Adonis / PERMANOVA**：檢定群間差異是否顯著。
* **ANOSIM / MRPP**：用於非參數群組差異檢定。

### 可視化方式：

* **2D/3D ordination 圖**：以點表示樣本，顏色代表組別。
* **群聚樹（Clustering tree）**：UPGMA 等法分群樣本。

### 目的：

* 分析 treatment/control 或不同環境來源間的群落組成差異。
* 評估樣本間是否存在明顯的菌群聚類趨勢。

---

## 📦 分析工具推薦

* `phyloseq`（R）：整合多樣性分析、可視化。
* `vegan`（R）：PERMANOVA、NMDS、距離矩陣分析。
* `qiime2 diversity core-metrics-phylogenetic`：提供一鍵分析。

---

這些多樣性分析能協助你解釋微生物群落的組成穩定性、環境變化影響、以及潛在的生物學意義。

若已完成 ASV 表（`seqtab.nochim`）與分類註解（`taxa`），即可導入 `phyloseq` 或 QIIME2 進行這些分析。
