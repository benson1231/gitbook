## Nextflow `params` 使用說明

在 Nextflow workflow 中，`params` 用來定義與接收使用者傳入的參數。
這些參數可以在 workflow 腳本中使用，或在執行時由命令列傳入。

---

### 🔹 定義 `params`

在 workflow 腳本（如 `main.nf`）中定義參數：

```groovy
params.input = "data/sample.csv"
params.outdir = "results"
```

可搭配 `params.<key>` 在流程中呼叫：

```groovy
process example {
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

### 🔹 使用命令列覆蓋 `params`

執行 workflow 時可用 `--<param>=<value>` 傳入參數：

```bash
nextflow run main.nf --input new_input.csv --outdir output_dir
```

---

### 🔹 結合 `nextflow.config`

可在 `nextflow.config` 檔案中預設參數：

```groovy
params {
  input = "data/default.csv"
  outdir = "results"
  threads = 4
}
```

這些預設值可以在命令列被覆蓋。

---

### 🔹 檢查參數值

可在流程前檢查參數內容：

```groovy
echo "Input file: ${params.input}"
echo "Output dir: ${params.outdir}"
```

---

`params` 是 Nextflow workflow 彈性化與模組化的重要工具，適用於設定輸入檔案、輸出資料夾、執行參數等。
