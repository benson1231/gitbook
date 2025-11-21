# Nextflow Channels 教學

Channels 是 Nextflow 的核心資料結構，用來實現 **資料驅動 (dataflow)** 的程式設計範式。它們用來在 processes 之間傳遞資料，或進行函數式的資料轉換。

---

## 1. Channel 類型

### Queue channel

* **FIFO**：先進先出 (First In, First Out)。
* **單向**：資料從生產者流向消費者。
* **非同步**：操作不會阻塞流程。
* 常見於 process output 或 `Channel.of()`、`Channel.fromPath()` 等 factory 建立。

```groovy
ch = Channel.of(1, 2, 3)
ch.view()
```

輸出：

```
1
2
3
```

### Value channel

* **單一值 (singleton)**。
* 可以無限次被讀取而不會消耗內容。
* 使用 `Channel.value()` 或 operators (`first`, `last`, `count`, `sum` 等) 建立。

```groovy
ch1 = Channel.of(1, 2, 3)
ch2 = Channel.value(1)

process SUM {
  input:
    val x
    val y
  output:
    stdout
  script:
    """
    echo $((x+y))
    """
}

workflow {
  SUM(ch1, ch2).view()
}
```

輸出：

```
2
3
4
```

---

## 2. Channel factories

### `value()`

建立 value channel。

```groovy
ch1 = Channel.value('Hello')
ch2 = Channel.value([1,2,3])
```

### `of()`

建立 queue channel。

```groovy
Channel.of(1, 3, 5, 7).view()
```

輸出：

```
1
3
5
7
```

> **註：** `Channel.from()` 已被棄用，建議使用 `Channel.of()`。

### `fromList()`

由 List 建立 channel。

```groovy
list = ['hello', 'world']
Channel.fromList(list).view()
```

輸出：

```
hello
world
```

### `fromPath()`

依檔案路徑建立 channel。

```groovy
Channel.fromPath('./data/meta/*.csv').view()
```

選項：

* `glob`: 是否使用通配符 (預設 true)
* `type`: file / dir / any
* `hidden`: 是否包含隱藏檔案
* `maxDepth`: 掃描子目錄層數
* `checkIfExists`: 檔案不存在是否拋錯

### `fromFilePairs()`

成對讀取檔案，輸出 tuple(key, \[file1, file2])。

```groovy
Channel.fromFilePairs('./data/ggal/*_{1,2}.fq').view()
```

輸出：

```
[liver, [liver_1.fq, liver_2.fq]]
[gut,   [gut_1.fq, gut_2.fq]]
```

### `fromSRA()`

直接從 **NCBI SRA** 抓取資料。

需要 `NCBI API key`：

```groovy
params.ncbi_api_key = '<Your API key>'

Channel.fromSRA(['SRP073307'], apiKey: params.ncbi_api_key).view()
```

輸出：

```
[SRR3383346, [SRR3383346_1.fastq.gz, SRR3383346_2.fastq.gz]]
...
```

---

## 3. 總結

* **Queue channel**：一次性 FIFO 資料流。
* **Value channel**：單一值，可重複讀取。
* **Channel factories**：

  * `value()`：建立單值 channel
  * `of()`：建立 FIFO channel（替代已棄用的 `from()`）
  * `fromList()`：由 List 建立
  * `fromPath()`：由檔案路徑建立
  * `fromFilePairs()`：配對檔案
  * `fromSRA()`：直接讀取 NCBI SRA

---

👉 建議搭配 Operators（如 `map`、`filter`、`splitText` 等）靈活使用，可大幅簡化 workflow 的資料處理。
