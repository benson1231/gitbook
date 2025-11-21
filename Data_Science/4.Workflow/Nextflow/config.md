## Nextflow `nextflow.config` 設定檔說明

`nextflow.config` 是 Nextflow workflow 的核心設定檔，用來管理全域參數、流程資源限制、執行後端（executor）、container 選項等。良好的 config 設計能提升 pipeline 的重現性、彈性與可維護性。

---

### 🔹 基本結構與語法

```groovy
params {
  input = 'data/samples.csv'
  outdir = 'results'
  threads = 4
}

process {
  cpus = 2
  memory = '4 GB'
  time = '2h'
}

docker {
  enabled = true
}
```

* `params`：定義 pipeline 參數，供 workflow 或 module 使用。
* `process`：定義流程的預設資源限制。
* `docker`/`singularity`：指定 container 執行設定。

---

### 🔹 Executor 執行後端選項

Nextflow 支援多種執行環境（executor）：

| Executor              | 描述                            |
| --------------------- | ----------------------------- |
| `local`               | 預設，使用本機執行流程                   |
| `slurm`               | 適用於 SLURM HPC 系統              |
| `pbs` / `sge`         | 適用於 PBS / Sun Grid Engine     |
| `lsf`                 | IBM LSF 系統                    |
| `ignite`              | 分散式叢集執行（實驗性）                  |
| `awsbatch`            | AWS Batch 雲端執行                |
| `google-lifesciences` | Google Cloud Life Sciences 執行 |
| `k8s`                 | Kubernetes 容器化環境              |

```groovy
process.executor = 'slurm'
process.queue = 'batch'
```

* 可用 `-profile` 自訂不同環境設定（見下）。

---

### 🔹 多重 profile 切換

可依執行環境建立多組設定，用 `-profile` 指定。

```groovy
profiles {
  local {
    process.executor = 'local'
  }
  docker {
    docker.enabled = true
    process.executor = 'local'
  }
  slurm {
    process.executor = 'slurm'
    singularity.enabled = true
  }
}
```

執行時切換：

```bash
nextflow run main.nf -profile slurm
```

---

### 🔹 模組化設定與外部 include

可透過 `includeConfig` 載入外部設定檔：

```groovy
includeConfig 'conf/genome.config'
```

* 適合將大型設定拆分成不同檔案。

---

### 🔹 使用 YAML 作為參數輸入

1. 建立 `params.yaml`：

```yaml
input: "data/samples.csv"
outdir: "results"
threads: 4
genome: "GRCh38"
```

2. 執行 workflow 並指定檔案：

```bash
nextflow run main.nf -params-file params.yaml
```

3. 在 workflow 中直接使用：

```groovy
params.input     // => "data/samples.csv"
params.genome    // => "GRCh38"
```

---

### 🔹 指定流程資源設定與選擇器（selector）

可針對個別流程名稱自訂資源：

```groovy
process {
  withName: fastqc {
    cpus = 4
    memory = '8 GB'
  }
  withName: 'bwa*' {
    cpus = 8
    time = '6h'
  }
  withName: /.*trim.*/ {
    cpus = 2
    memory = '3 GB'
  }
}
```

* 支援文字、萬用字元、正規表示式作為條件。

---

### 🔹 動態資源分配（Dynamic assignment）

流程內可使用 `params` 或 `task.attempt` 進行動態設定：

```groovy
process ALIGN {
  cpus = params.threads
  memory = { 2.GB * task.attempt }
  time = { task.attempt <= 2 ? '2h' : '4h' }
}
```

* 可根據執行次數動態增加資源。

---

### 🔹 錯誤策略與容錯控制

可設定流程的錯誤處理行為：

```groovy
process {
  errorStrategy = 'retry'
  maxRetries = 3
  maxErrors = 10
}
```

* `errorStrategy` 可為：`retry`, `finish`, `ignore`, `terminate`
* `maxRetries`：流程最多可重跑次數。
* `maxErrors`：所有流程最大錯誤數。

---

### 🔹 實用功能與補充

* `queueSize`：同時執行的流程數量限制（對 executor 有效）
* `scratch`：指定暫存目錄（部分 cluster 支援）
* `cache`：控制 process cache 行為（開關、路徑）
* `singularity.enabled = true`：Singularity container 支援
* `wave.enabled = true`：支援 Wave container service
* `dag.file = 'dag.png'`：輸出流程圖

---

良好撰寫的 `nextflow.config` 能使你的 pipeline 在本機、HPC 或雲端間無縫切換，並方便團隊協作與部署。建議搭配 profile 與 YAML 管理參數以強化模組化與重現性。
