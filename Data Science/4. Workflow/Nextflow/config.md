## Nextflow `nextflow.config` 設定檔說明

`nextflow.config` 是 Nextflow pipeline 的全域設定檔，用來指定參數、執行環境、container、流程資源限制等。

---

### 🔹 基本結構

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

---

### 🔹 executor 支援的選項

Nextflow 支援多種執行後端（executor），可根據執行環境選擇不同模式：

| Executor              | 描述                     |
| --------------------- | ---------------------- |
| `local`               | 預設本地執行，使用單一機器的 CPU/記憶體 |
| `slurm`               | 適用於 HPC 系統（如 SLURM 集群） |
| `pbs`                 | 適用於 PBS/Torque 系統      |
| `sge`                 | 適用於 Sun Grid Engine    |
| `lsf`                 | IBM Platform LSF       |
| `ignite`              | 分散式叢集執行（experimental）  |
| `awsbatch`            | 在 AWS Batch 上執行流程      |
| `google-lifesciences` | 在 Google Cloud 上執行     |
| `k8s`                 | 在 Kubernetes 上執行       |
| `process.executor`    | 可在 profile 或流程中個別指定    |

**範例（設定 SLURM）：**

```groovy
process.executor = 'slurm'
process.queue = 'batch'
```

你也可以在 profile 區塊中切換不同 executor 設定（見下方）。

---

### 🔹 多 profile 設定

可透過 `-profile` 切換不同環境或執行方式：

```groovy
profiles {
  standard {
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

執行時指定 profile：

```bash
nextflow run main.nf -profile docker
```

---

### 🔹 module 與 workflow 預設值設定（DSL2）

```groovy
params {
  genome = 'GRCh38'
}

includeConfig 'conf/genome.config'
```

---

### 🔹 使用外部 YAML 檔讀取參數

Nextflow 支援從 YAML 檔案讀取參數並傳入 workflow：

1. 建立一個 `params.yaml`：

```yaml
input: "data/samples.csv"
outdir: "results"
threads: 4
genome: "GRCh38"
```

2. 執行 workflow 並載入 YAML：

```bash
nextflow run main.nf -params-file params.yaml
```

3. 在 workflow 中使用這些參數：

```groovy
params.input     // -> "data/samples.csv"
params.genome    // -> "GRCh38"
```

---

### 🔹 自訂參數與流程共用

config 中設定的 `params` 可在任何流程或 module 中使用：

```groovy
input:
  val genome from params.genome
```

---

### 🔹 小技巧與額外功能

* `executor.cpus`, `executor.memory`, `executor.queue` 可依環境自定。
* `withName: <process_name>` 可指定個別流程資源：

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

* 支援 selector 表達式：

  * 字串：`withName: 'process_name'`
  * 萬用字元：`withName: 'align*'`
  * 正規表示式：`withName: /.*qc$/`

---

### 🔹 動態資源設定（dynamic resource assignment）

你可以在 `process` 中使用變數（如 `params.threads`）來動態設定資源：

```groovy
process ALIGN {
  cpus = params.threads
  memory = { 2.GB * task.attempt }
  time = { task.attempt <= 2 ? '2h' : '4h' }
}
```

* `params.threads`：來自 config 或 YAML 的外部參數
* `task.attempt`：表示第幾次重試，可用於設定遞增的時間或記憶體
* 動態計算可用 `{}` 包住表達式

---

### 🔹 資源限制與保護（resource limits）

可透過 `errorStrategy` 與 `maxRetries`, `maxErrors` 控制流程容錯：

```groovy
process {
  errorStrategy = 'retry'
  maxRetries = 3
  maxErrors = 10
}
```

* `retry`：遇錯重跑；也可設為 `finish`, `ignore`, `terminate`
* `maxRetries`：最多重試幾次
* `maxErrors`：流程允許總錯誤數上限（超過會終止）

---

良好地撰寫 `nextflow.config` 並結合外部 YAML 管理參數，有助於讓 workflow 更可重複、跨平台適用與團隊共享。
