## 如何執行 Nextflow 與常用參數說明

Nextflow 是一個針對生物資訊與資料分析工作流程設計的工作流程語言，透過以下步驟與參數說明，可以有效執行與掌控流程行為。

---

### ✅ 基本執行指令

```bash
nextflow run <workflow.nf> [options]
```
或
```bash
nextflow run <github-repo> [-r <revision>] [options]
```

範例：
```bash
nextflow run main.nf -profile docker --input data/sample.csv --outdir results
```

---

### 🔹 常用參數總覽

| 類型       | 參數                           | 說明 |
|------------|--------------------------------|------|
| Pipeline   | `<script.nf>` / `<repo>`       | 本地 `.nf` 檔案或 GitHub repo，如 `nf-core/rnaseq` |
| Profile    | `-profile <name>`              | 使用不同的配置檔，如 `docker`, `singularity`, `slurm` |
| 修復/重跑  | `-resume`                       | 繼續先前未完成流程，使用 cache |
| 修訂版     | `-r <version>`                  | 指定 pipeline 的 Git 版本或 tag，例如 `-r 3.9` |
| Config     | `-c <config file>`              | 指定額外 config 設定檔 |
| Params     | `--param <value>`               | 自訂流程參數，例如 `--genome GRCh38` |
| Params 檔 | `-params-file <yaml/json>`      | 使用 YAML/JSON 格式傳入大量參數 |
| 輸出路徑   | `--outdir <path>`               | 指定輸出結果資料夾 |
| 後台執行   | `-bg`                           | 背景執行流程 |
| Dry run    | `-n` 或 `-dry-run`              | 模擬流程，不實際執行命令（可搭配 `-preview`） |
| 工作路徑   | `-work-dir <path>`              | 自訂中間檔存放目錄（預設為 `./work`） |
| 顯示控制   | `-ansi-log false`               | 關閉彩色 log（適用於無支援 ANSI 終端） |

---

### 🔹 常見執行情境

#### ✅ 執行本地流程（使用 docker profile）：
```bash
nextflow run main.nf -profile docker --input samples.csv --fasta hg38.fa --outdir results
```

#### ✅ 執行 nf-core 流程（使用 Singularity 且指定版本）：
```bash
nextflow run nf-core/rnaseq -profile singularity -r 3.13 --input samplesheet.csv --genome GRCh38
```

#### ✅ 使用 YAML 傳入參數：
```bash
nextflow run main.nf -profile docker -params-file params.yaml
```
```yaml
# params.yaml 範例
input: "data/sheet.csv"
outdir: "results"
genome: "GRCh37"
```

#### ✅ 指定額外 config（如 SLURM 設定）：
```bash
nextflow run main.nf -c conf/slurm.config -profile slurm
```

#### ✅ 重跑未完成流程：
```bash
nextflow run main.nf -resume
```

---

### 🧪 查看流程版本與資源記錄

```bash
nextflow log              # 查看所有歷史流程紀錄
nextflow log <run-name>  # 查看特定流程細節
```

---

### 📌 小技巧

- 可用 `nextflow config` 查看當前 config 設定細節。
- 使用 `NXF_OPTS` 限制記憶體：
```bash
export NXF_OPTS="-Xms1g -Xmx4g"
```
- 配合 `screen` 或 `tmux` 可避免中斷。
- 常與 GitHub Actions 或 SLURM 結合做 CI/CD 與平行分析。

---

如需更多範例，可參考官方文件：[https://www.nextflow.io/docs/latest/](https://www.nextflow.io/docs/latest/) 。
