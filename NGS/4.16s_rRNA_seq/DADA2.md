# DADA2

**DADA2 (Divisive Amplicon Denoising Algorithm 2)** 是一套專為 16S rRNA 與其他擴增子定序設計的錯誤修正與序列推斷工具。其核心理念為：

> 直接建模定序錯誤機率，從原始 reads 中推斷真實存在的變異序列（ASVs, Amplicon Sequence Variants），而非依靠傳統聚類產生 OTUs。

---

## 🔬 核心流程

1. **品質過濾與裁剪 (filterAndTrim)**  
   - 移除低品質 reads，截斷低品質尾端，去除非生物學序列 (primers/adapters)。

2. **學習錯誤模型 (learnErrors)**  
   - 建立鹼基錯誤率矩陣，捕捉平台特異的錯誤模式。

3. **去噪 (dada)**  
   - 應用錯誤模型，區分真實變異與錯誤產生的假變異。

4. **配對合併 (mergePairs)**  
   - 將正反向 reads 對齊合併，保留一致序列，提升覆蓋長度與準確度。

5. **去除嵌合體 (removeBimeraDenovo)**  
   - 移除 PCR 擴增過程形成的假序列。

6. **建立 ASV 表 (makeSequenceTable)**  
   - 生成樣本 × ASV 的 abundance 矩陣。

7. **分類註解 (assignTaxonomy / addSpecies)**  
   - 透過參考資料庫 (如 SILVA、RDP) 對 ASVs 進行分類標註。

---

## ✅ DADA2 相較 OTU 的優勢

| 面向       | OTU 方法 (97% 相似度) | DADA2 / ASV 方法 |
|------------|--------------------|------------------|
| **解析度**   | 近似群體，無法區分單一鹼基差異 | 單一鹼基解析度 (更精確) |
| **依賴參考** | 可依賴或不依賴           | 僅在註解階段依賴 |
| **重現性**   | 聚類參數影響大，不穩定       | 高度可重現 (ASV 為唯一單位) |
| **生物學解釋** | 低 (混合菌群)             | 高，可追蹤特定變異 |

---

## 🧬 輸出成果

- `seqtab.nochim`：ASV abundance 表格 (樣本 × ASV)
- `taxa`：ASV 對應的分類階層 (Kingdom → Genus/Species)

這些結果可直接匯入 **phyloseq** 進行：
- α 多樣性與 β 多樣性分析
- 統計檢定與可視化 (richness, ordination, barplot)

---

## 📖 延伸閱讀
- [DADA2 官方文件](https://benjjneb.github.io/dada2/index.html)
- [MiSeq SOP 教學](https://mothur.org/wiki/miseq_sop/)
- [phyloseq 套件](https://joey711.github.io/phyloseq/)
