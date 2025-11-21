# vcf2maf 使用速查表

**更新日期**: 2025.8.22
**參考來源**: [https://github.com/mskcc/vcf2maf](https://github.com/mskcc/vcf2maf)

---

## 前言

**vcf2maf** 是由 MSKCC（Memorial Sloan Kettering Cancer Center）開發的工具，用於將 VCF（Variant Call Format）轉換為 MAF（Mutation Annotation Format）。

在癌症基因體學研究中，VCF 是變異檔案的標準格式，但其結構較複雜、不利於下游統計分析與可視化。而 MAF 則是專為癌症突變資料設計的表格式檔案，廣泛應用於 TCGA 與多數癌症基因體研究。**vcf2maf** 結合 **Ensembl VEP (Variant Effect Predictor)** 進行註解，將 VCF 轉換為標準化的 MAF，方便後續的突變分析與統計工作。

主要用途包括：

* 將 VCF 格式轉換為 MAF 格式
* 利用 VEP 註解變異的功能影響與基因資訊
* 與下游癌症突變分析流程（如 maftools）相整合

---

## 一、基本使用範例

```bash
# 將 VCF 轉換為 MAF
perl vcf2maf.pl \
  --input-vcf example.vcf \
  --output-maf example.maf \
  --vep-path /path/to/vep \
  --vep-data /path/to/vep_data

# 指定腫瘤與正常樣本名稱
perl vcf2maf.pl \
  --input-vcf example.vcf \
  --output-maf example.maf \
  --tumor-id TumorSample1 \
  --normal-id NormalSample1
```

---

## 二、常用參數

* **`--input-vcf`**：輸入 VCF 檔案。
* **`--output-maf`**：輸出 MAF 檔案。
* **`--tumor-id`**：指定腫瘤樣本名稱（必填）。
* **`--normal-id`**：指定正常樣本名稱（可選）。
* **`--ref-fasta`**：參考基因組 FASTA 檔案（與 VCF 一致）。
* **`--vep-path`**：VEP 程式路徑。
* **`--vep-data`**：VEP 註解資料目錄。
* **`--ncbi-build`**：基因組版本（例如 GRCh37 或 GRCh38）。
* **`--retain-info`**：指定保留 VCF 中的 INFO 欄位。
* **`--filter-vcf`**：過濾 VCF 中不符合條件的變異。

---

## 三、輸出檔案

生成的 **MAF 檔案** 為制式化表格，包含以下資訊：

* 基因名稱（Hugo\_Symbol）
* 染色體、位置、REF/ALT
* 變異型別（SNV、INDEL）
* 功能註解（Variant\_Classification, Variant\_Type）
* 腫瘤樣本、正常樣本的基因型資訊

範例：

```
Hugo_Symbol  Chromosome  Start_Position  End_Position  Variant_Classification  Tumor_Sample_Barcode
TP53         17          7579472         7579472       Missense_Mutation       TumorSample1
```

---

## 四、常見應用

* 整合 TCGA 或 ICGC 的突變資料進行比較研究
* 搭配 **maftools** 進行可視化與統計分析
* 作為腫瘤基因組 pipelines 的標準輸出格式

---

## 五、注意事項

* **需安裝 VEP**，並確保版本與 vcf2maf 相容。
* VCF 的參考基因組（ref fasta）需與 VEP 註解資料一致。
* MAF 檔案體積可能較大，建議壓縮保存。
* vcf2maf 是針對癌症研究而設計，若用於其他研究需注意格式適配性。

---

**vcf2maf** 是癌症基因體學研究中，將 VCF 統一轉換與註解為 MAF 格式的標準工具，能大幅簡化資料前處理並支援後續多種分析流程。
