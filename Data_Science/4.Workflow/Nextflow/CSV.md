## 使用 Channel 處理 CSV 檔案

在 Nextflow 中，可以通過 `Channel.fromPath` 搭配 `splitCsv()` 與 `map` 方法來讀取與處理 `.csv` 檔案內容。

---

### 🔹 讀取 CSV 檔案

使用 `fromPath` 搭配 `splitCsv()` 將 CSV 拆分為結構化資料：

```groovy
Channel.fromPath('samplesheet.csv')
       .splitCsv(header: true)
       .view()
```

* `header: true` 表示第一列為欄位名稱。
* 輸出為 `Map`，欄位名稱作為 key。

---

### 🔹 搭配 `map` 處理資料

可對每一筆 row 資料進行轉換或提取特定欄位：

```groovy
Channel.fromPath('samplesheet.csv')
       .splitCsv(header: true)
       .map { row -> tuple(row.sample_id, file(row.fastq)) }
       .set { sample_ch }
```

* `row` 為每筆資料（Map 結構），例如 `row.sample_id`、`row.fastq`。
* `tuple(...)` 用於輸出成對資料給 process 使用。

---

### 🔹 於 Process 中使用

```groovy
process qc {
  input:
  tuple val(sample_id), path(reads)

  script:
  """
  echo Running QC for $sample_id with $reads
  """
}
```

---

此方法常用於從 metadata 表或樣本表中讀取參數，方便 batch 處理與流程動態化設計。