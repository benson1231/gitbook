## Nextflow Module（模組）設計教學

Nextflow DSL2 支援 module 化流程設計，使 pipeline 更具重用性與可維護性。模組通常定義於 `modules/` 資料夾中，並可在主 workflow 中呼叫。

---

### 🔹 基本結構

一個 module 是一個獨立的 `process`，可包含 `input:`, `output:`, `script:` 等區塊。

**範例：`modules/qc/fastqc.nf`**

```groovy
process FASTQC {
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)

  output:
    path "${sample_id}_fastqc.zip" emit: fastqc_zip

  script:
  """
  fastqc $reads -o .
  mv *.zip ${sample_id}_fastqc.zip
  """
}
```

---

### 🔹 在主 workflow 中使用 module

**主腳本：`main.nf`**

```groovy
workflow {
  Channel.fromPath('data/*.fq.gz')
         .map { file -> tuple(file.baseName, file) }
         .set { read_ch }

  FASTQC(read_ch)
  FASTQC.out.fastqc_zip.view()
}
```

---

### 🔹 `modules.json` 或 `modules/` 管理

你可以用以下方式引用模組：

```groovy
include { FASTQC } from './modules/qc/fastqc.nf'
```

或使用 DSL2 的 modules 匯入機制：

```groovy
include { FASTQC } from './modules/qc'         // 需模組名稱與檔名一致
```

---

### 🔹 模組設計原則

* 每個模組應保持功能單一
* 使用 `emit:` 命名輸出，方便 workflow 使用
* 可與 `nextflow config` 結合設定預設參數

---

模組化可提升 Nextflow 專案的可讀性與可維護性，特別適用於大型 pipeline 或多人開發環境。
