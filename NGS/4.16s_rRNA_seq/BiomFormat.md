# 🧬 BIOM 格式介紹

## 📌 什麼是 BIOM 格式？

BIOM (**Biological Observation Matrix**) 是一種專門為微生物社群分析設計的檔案格式。它主要用來儲存 OTU（Operational Taxonomic Unit）或 ASV（Amplicon Sequence Variant）在不同樣本中的豐度資訊，以及相關的 **樣本與分類註解**。

BIOM 格式最初由 QIIME 社群提出，後來成為微生物學與生態學研究中的標準格式。

---

## 📂 BIOM 的核心特徵

* **稀疏矩陣結構**：適合存放龐大且稀疏的 OTU/ASV table。
* **支援多種 metadata**：可同時包含樣本 metadata（如取樣時間、組織來源）、OTU/ASV metadata（如 taxonomy）。
* **跨平台兼容**：廣泛支援於 QIIME、Mothur、PICRUSt、phyloseq 等工具。
* **JSON / HDF5 格式**：

  * BIOM v1 → 基於 JSON（文字格式，適合小型資料）
  * BIOM v2 → 基於 HDF5（二進位格式，適合大型資料）

---

## 📊 BIOM 格式結構

BIOM table 基本上是一個「樣本 × OTU/ASV」矩陣，外加 metadata，例如：

```
# Constructed from biom file
#OTU ID   Sample1   Sample2   Sample3   taxonomy
OTU_1     10        0         2         k__Bacteria; p__Firmicutes; g__Lactobacillus
OTU_2     0         5         1         k__Bacteria; p__Actinobacteria; g__Bifidobacterium
```

這樣的 table 會在 BIOM 檔案中以壓縮的 JSON/HDF5 儲存，並能保存更豐富的 metadata。

---

## 🔧 常見用途

1. **QIIME / QIIME2**：輸入輸出 OTU/ASV 表格的標準格式。
2. **PICRUSt / PICRUSt2**：必須使用 BIOM 格式作為輸入。
3. **phyloseq (R package)**：可直接匯入 BIOM 格式建立 phyloseq object。

---

## 🔄 BIOM 與文字格式的轉換

安裝 `biom-format` 工具後，可自由轉換：

```bash
# 將 tab-delimited 文字檔轉為 BIOM 格式
biom convert -i otu_table.txt -o otu_table.biom --table-type="OTU table" --to-hdf5

# 將 BIOM 檔轉回文字檔
biom convert -i otu_table.biom -o otu_table.txt --to-tsv --header-key taxonomy
```

---

## 📌 總結

* BIOM 格式是微生物分析的核心資料結構，用來存放 **豐度表 + 註解**。
* 它的優勢在於 **壓縮、支援 metadata、跨平台使用**。
* 幾乎所有常見的 16S/微生物分析工具（QIIME, PICRUSt, phyloseq）都能處理 BIOM 格式。
