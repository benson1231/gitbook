## Nextflow 常用指令介紹

以下是 Nextflow 中 `run`、`log` 與 `clean` 指令的簡要說明與用法。

---

### 🔹 `nextflow run`

執行 Nextflow workflow 的主要指令。

**語法：**

```bash
nextflow run <workflow_path_or_repo> [options]
```

**常用選項：**

* `-profile <name>`：使用指定的設定檔（如 `docker`, `conda`, `slurm` 等）。
* `--<param>=<value>`：傳遞參數給 workflow，例如 `--input samplesheet.csv`。
* `-resume`：從先前失敗的地方繼續執行。
* `-name <name>`：為此次執行命名，便於追蹤。
* `-ansi-log false`:關閉 Nextflow 預設的彩色日誌輸出（ANSI 格式）

**範例：**

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

---

### 🔹 `nextflow log`

顯示 workflow 的執行記錄。

**語法：**

```bash
nextflow log [options]
```

**常用選項：**

* `-f <fields>`：指定要輸出的欄位，例如 `status,name,duration`。
* `-last`：顯示最近一次的執行紀錄。

**範例：**

```bash
nextflow log -last -f 'runName,status,duration'
```

---

### 🔹 `nextflow clean`

移除 workflow 執行中產生的暫存資料。

**語法：**

```bash
nextflow clean [options]
```

**常用選項：**

* `-n` 或 `--dry-run`：僅顯示會被刪除的檔案，不實際執行刪除。
* `-f` 或 `--force`：強制執行刪除，不再提示確認。

**範例：**

```bash
nextflow clean -f     # 強制清除所有暫存
nextflow clean -n     # 預覽將刪除的內容
```

---

這些指令有助於有效管理 Nextflow workflow 的執行與資源清理，可搭配 `nextflow config`、`nextflow info` 等其他指令進行更完整的流程控制與除錯。
