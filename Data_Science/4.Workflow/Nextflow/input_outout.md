## Nextflow Process I/O 與控制屬性說明

在 Nextflow 的 `process` 區塊中，`input:` 與 `output:` 是資料傳遞的核心，而額外屬性如 `script`、`cpus`、`errorStrategy` 等則定義流程執行細節。熟練使用這些設定能有效提升流程模組化與平行效率。

---

### 🔹 `input:` 區塊（輸入）

常見輸入型態：

| 類型      | 說明                        |
| ------- | ------------------------- |
| `val`   | 純量變數（如 string、int）        |
| `path`  | 檔案或資料夾（會自動 staged 至工作目錄）  |
| `tuple` | 多項輸入組合                    |
| `each`  | 展開 list 為多個 item 傳入（平行處理） |
| `stdin` | 來自標準輸入的文字流                |

**範例：**

```groovy
input:
  val sample_id
  path "${sample_id}.fq"
```

```groovy
input:
  tuple val(id), path(fq)
```

```groovy
input:
  each val(name) from Channel.of('A','B','C')
```

---

### 🔹 `output:` 區塊（輸出）

常見輸出型態：

| 類型       | 說明                |
| -------- | ----------------- |
| `path`   | 匹配的輸出檔案           |
| `val`    | 單一變數輸出（如純量結果）     |
| `tuple`  | 同時輸出多項元素          |
| `emit:`  | 命名輸出，用於 DSL2 模組引用 |
| `stdout` | 將 stdout 作為輸出     |
| `stderr` | 將 stderr 作為輸出     |
| `into:`  | 將輸出導入其他 channel   |

**範例：**

```groovy
output:
  path "*.bam" into bam_ch
```

```groovy
output:
  tuple val(sample_id), path("*.bam") emit: aligned
```

```groovy
output:
  path "*.log" collect: true into log_ch
```

---

### 🔹 輸出收集與包裝

* `collect: true`：將多個檔案合併為一個 list 輸出
* `flatten: true`：攤平成單一串流（去除巢狀）

```groovy
output:
  path "*.txt" collect: true emit: reports
```

---

### 🔹 多輸入設計與組合

可使用 `combine()` 或 `join()` 將不同來源 channel 整合後輸入：

```groovy
ids.combine(reads).set { input_ch }
```

```groovy
input:
  tuple val(id), path(fq) from input_ch
```

---

### 🔹 額外流程屬性

| 屬性               | 說明                                           |
| ---------------- | -------------------------------------------- |
| `script:`        | 實際要執行的命令或 bash script                        |
| `cpus`, `memory` | 指定流程所需資源，支援動態計算（如 `{ 2.GB * task.attempt }`） |
| `time`           | 最長執行時間限制，如 `'2h'`, `'30min'`                 |
| `container`      | 指定特定 container image                         |
| `publishDir`     | 輸出檔案自動複製到指定資料夾                               |
| `label`          | 可與 config selector 搭配調控資源                    |
| `tag`            | 加入標籤供 log 與 UI 顯示                            |
| `errorStrategy`  | 遇錯行為，如 `'retry'`, `'finish'`, `'ignore'`     |
| `maxRetries`     | 最多允許重跑次數                                     |
| `when:`          | 條件式執行（true 才會執行此 process）                    |

---

### 🔹 簡易範例：整合所有設定

```groovy
process ALIGN {
  tag "$sample_id"
  label 'bwa_align'

  input:
    tuple val(sample_id), path(reads)

  output:
    path "${sample_id}.bam" emit: aligned

  script:
  """
  bwa mem ref.fa $reads > ${sample_id}.bam
  """
}
```

---

良好的 I/O 與流程屬性設計，是建立彈性、高效、可維護 workflow 的關鍵。建議搭配 `params`、`channel` 操作與 config selector 強化模組整合性。
