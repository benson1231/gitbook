# SAM Format

SAM（Sequence Alignment/Map）格式是次世代定序（NGS）中記錄序列對齊結果的純文字格式，適用於記錄將定序讀段（read）對齊至參考基因組後的資訊，是 BAM 格式的原始形式。

---

## 一、SAM 格式基本結構

SAM 檔案由兩部分組成：

1. **Header（標頭）**：以 `@` 開頭，描述參考基因組與對齊參數等資訊。例如：

   ```
   @SQ	SN:chr1	LN:248956422
   @PG	ID:bwa	PN:bwa	VN:0.7.17
   ```

2. **Alignment（對齊紀錄）**：每行為一筆讀段對齊結果，共包含至少 11 個欄位：

   | 欄位    | 說明                        |
   | ----- | ------------------------- |
   | QNAME | 讀段名稱（query template name） |
   | FLAG  | 旗標，描述讀段的對齊狀態與方向等          |
   | RNAME | 對齊之參考序列名稱（如 chr1）         |
   | POS   | 對齊起始位置（從 1 開始）            |
   | MAPQ  | 對齊品質分數（mapping quality）   |
   | CIGAR | 對齊摘要（表示插入、缺失、比對長度等）       |
   | RNEXT | 配對讀段的參考名稱（若無配對，為 \*）      |
   | PNEXT | 配對讀段的起始位置                 |
   | TLEN  | 模板長度（template length）     |
   | SEQ   | 讀段序列                      |
   | QUAL  | 品質分數（Phred 品質，以 ASCII 表示） |

   **註**：之後可接額外的 TAG 欄位（以 `TAG:TYPE:VALUE` 格式表示）

---

## 二、FLAG 欄位說明

FLAG 為整數數值，經位元運算表示多種讀段特性，可用工具如 `samtools flagstat` 解析。

常見 FLAG 內容如下：

| FLAG 值 | 意義                       |
| ------ | ------------------------ |
| 1      | 讀段為配對（paired-end）        |
| 2      | 對齊成功且成對（properly paired） |
| 4      | 未對齊                      |
| 16     | 讀段在負向鏈對齊                 |
| 32     | 配對讀段在負向鏈                 |
| 64     | 此為第一讀段（first in pair）    |
| 128    | 此為第二讀段（second in pair）   |

---

## 三、SAM 格式應用與工具

* **檢查比對品質**：`samtools view`, `samtools flagstat`
* **格式轉換**：轉為 BAM 使用 `samtools view -b`
* **視覺化對齊結果**：轉為 BAM 並建立索引後載入 IGV

---

## 四、注意事項

* SAM 為純文字格式，檔案體積較大，分析建議轉為 BAM 處理
* 檔案內容可用文字編輯器開啟檢查錯誤或格式問題

---

SAM 是對齊結果的標準格式，易於閱讀與除錯，適合開發與測試階段使用，最終則多轉為 BAM 以利儲存與分析。
