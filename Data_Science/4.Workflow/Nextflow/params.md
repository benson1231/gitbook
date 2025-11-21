## Nextflow `params` 使用說明

在 Nextflow workflow 中，`params` 是用來定義與接收使用者傳入的參數機制。這些參數可以在 workflow 腳本中設定預設值，或在執行時透過命令列動態指定。

---

### 🔹 在 `main.nf` 中定義參數

```groovy
params.input = "data/sample.csv"
params.outdir = "results"
params.threads = 4
```

在流程中可透過 `params.<key>` 使用：

```groovy
process EXAMPLE {
  input:
    path sample_file from file(params.input)

  output:
    path "*.txt" into result_files

  script:
  """
  cp $sample_file output.txt
  """
}
```

---

### 🔹 執行時使用命令列覆蓋 `params`

```bash
nextflow run main.nf --input new_input.csv --outdir output_dir --threads 8
```

---

### 🔹 在 `nextflow.config` 中設定預設值

```groovy
params {
  input = "data/default.csv"
  outdir = "results"
  threads = 4
  name = "User"
}
```

這些值會被命令列參數覆蓋。

---

### 🔹 使用三元運算設定預設值

可依據條件決定要使用預設值或外部指定：

```groovy
// 若 params.name 有值，使用它；否則使用 'Nextflow'
def names_ch = params.name ? Channel.of(params.name) : Channel.of('Nextflow')
```

---

### 🔹 結合 YAML 檔傳入參數

Nextflow 支援從外部 YAML 載入參數：

📄 `params.yaml`

```yaml
input: "data/sample.csv"
outdir: "results"
threads: 8
```

執行方式：

```bash
nextflow run main.nf -params-file params.yaml
```

---

### 🔹 檢查與除錯 `params`

在 workflow 或流程中直接列印參數：

```groovy
echo "Input file: ${params.input}"
echo "Output dir: ${params.outdir}"
```

---

### 🔹 常見使用場景

| 場景         | 範例                                  |
| ---------- | ----------------------------------- |
| 指定輸入資料     | `--input data/sample1.csv`          |
| 設定輸出目錄     | `--outdir output/`                  |
| 控制執行緒數量    | `--threads 8`                       |
| 選擇分析基因組版本  | `--genome GRCh38`                   |
| 指定模式或行為控制值 | `--mode full`（可配合 `params.mode` 判斷） |

---

良好地使用 `params` 有助於讓 workflow 更具彈性與可重複性，便於跨平台部署與多人合作，亦是模組化與可參數化設計的基礎。
