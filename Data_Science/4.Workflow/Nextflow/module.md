## Nextflow Module（模組）設計教學

Nextflow DSL2 支援模組化流程設計，使 pipeline 更具重用性、可測試性與可維護性。模組通常定義於 `modules/` 或 `subworkflows/` 資料夾中，並在主 workflow 中引入。

---

### 🔹 基本模組結構（單一 process 模組）

一個模組通常是一個 `process` 定義，包含 `input:`、`output:`、`script:` 區塊，並使用 `emit:` 命名輸出。

📄 `modules/qc/fastqc.nf`

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

### 🔹 主 workflow 中呼叫模組

📄 `main.nf`

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

### 🔹 模組匯入方式（include）

在 `main.nf` 中使用 `include` 引入模組：

```groovy
include { FASTQC } from './modules/qc/fastqc.nf'
```

若模組名稱與檔名一致，也可只給路徑：

```groovy
include { FASTQC } from './modules/qc'
```

若有多個模組要引用，可批次包含：

```groovy
include { FASTQC; TRIMMOMATIC; BWA } from './modules/qc'
```

---

### 🔹 使用 emit 與輸出名稱

模組的輸出應透過 `emit:` 來命名，使主 workflow 可明確引用：

```groovy
output:
  path "${sample_id}.bam" emit: aligned_bam
```

```groovy
workflow {
  BWA(indexed_reads)
  BWA.out.aligned_bam.view()
}
```

---

### 🔹 模組化設計原則

| 原則     | 說明                                 |
| ------ | ---------------------------------- |
| 單一功能原則 | 每個模組只負責一個明確步驟（如 fastqc, trim, map） |
| 輸入輸出明確 | 使用 `tuple` 與 `emit` 定義清楚的輸入與輸出     |
| 可重用性   | 模組應盡量與特定資料路徑或參數脫鉤，可透過 `params` 設定  |
| 無副作用   | 避免模組修改共用檔案，所有輸出皆在工作目錄中產生           |

---

### 🔹 進階：模組參數與 config 整合

模組中可使用 `params` 提供額外參數支援：

```groovy
params {
  fastqc_threads = 2
}
```

模組中：

```groovy
script:
"""
fastqc --threads ${params.fastqc_threads} $reads -o .
"""
```

---

### 🔹 子流程模組（subworkflow）

除了單一 `process`，也可將多個模組組成 `subworkflow`，以便封裝高階邏輯：

📄 `subworkflows/qc_pipeline.nf`

```groovy
workflow QC_PIPELINE {
  take:
    reads

  main:
    TRIMMOMATIC(reads)
    FASTQC(TRIMMOMATIC.out.trimmed_reads)

  emit:
    trimmed = TRIMMOMATIC.out.trimmed_reads
    report  = FASTQC.out.fastqc_zip
}
```

主流程中：

```groovy
include { QC_PIPELINE } from './subworkflows/qc_pipeline.nf'
QC_PIPELINE(reads_ch)
```

---

模組化是大型 Nextflow workflow 的關鍵，配合 `emit:`、`include`、`params` 等機制可打造具可讀性、可測試性、可版本化的流程元件。建議搭配 `nf-core` 的模組撰寫標準與自動化測試，強化流程穩定性與共享性。
