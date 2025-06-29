## Nextflow I/O（輸入與輸出）說明

Nextflow 使用 `input:` 和 `output:` 區塊在 `process` 中定義資料的接收與傳出方式。這些區塊可以搭配多種資料型態與 channel 操作，支援高效且彈性的流程設計。

---

### 🔹 `input:` 區塊

常見輸入類型：

* `val`：傳入純量（string, int, 等）
* `path`：傳入檔案或目錄
* `tuple`：同時傳入多個值

**範例：**

```groovy
input:
  val sample_id
  path reads
```

或：

```groovy
input:
  tuple val(sample_id), path(reads)
```

---

### 🔹 `output:` 區塊

常見輸出類型：

* `path`：將產生的檔案輸出為 Channel
* `val`：輸出變數
* `tuple`：多輸出值的組合
* `emit:`：命名輸出，供 DSL2 module 使用

**範例：**

```groovy
output:
  path "*.bam" into bam_ch
```

或（多輸出）：

```groovy
output:
  tuple val(sample_id), path("*.bam") into bam_ch
```

---

### 🔹 多重輸出與 `collect`

若流程會產出多個檔案或變數並希望一次輸出：

**範例：**

```groovy
output:
  path "results/*.txt" collect: true into txt_ch
```

* `collect: true` 會將多個檔案打包為一個 list 輸出

---

### 🔹 多輸入流程設計

**組合多個 channel：**

```groovy
Channel.from(['sample1', 'sample2']).set { ids }
Channel.fromPath('*.fq').set { files }

ids.combine(files).set { input_ch }

process run {
  input:
  tuple val(id), path(fq)
  ...
}
```

---

Nextflow 的 I/O 機制設計靈活，透過 tuple、combine、collect 等操作，可達成模組化、平行化與高效率的流程建構。
