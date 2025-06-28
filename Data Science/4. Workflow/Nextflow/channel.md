## Nextflow `channel` 使用說明

在 Nextflow 中，`Channel` 是資料傳遞的核心機制。流程之間的資料流通都是透過 Channel 來實現的。

---

### 🔹 建立 Channel

**1. 從檔案建立：**

```groovy
Channel.fromPath('data/*.fastq')
```

**2. 從清單建立：**

```groovy
Channel.from(['sample1', 'sample2', 'sample3'])
```

**3. 使用 `Channel.of()` 明確建立：**

```groovy
Channel.of('a', 'b', 'c')
```

* `of` 是 `from` 的別名，效果相同但語義更簡潔。
* 適用於明確列出靜態元素。

**4. 空 Channel（供動態填入）：**

```groovy
Channel.empty()
```

---

### 🔹 Channel 操作

**1. `map`：轉換資料內容**

```groovy
ch_names = Channel.from(['file1', 'file2']).map { it.toUpperCase() }
```

**2. `filter`：篩選資料**

```groovy
ch_filtered = ch_files.filter { it.endsWith('.fastq') }
```

**3. `set` / `tuple`：建立複合資料流**

```groovy
Channel.of('sample1', 'sample2')
      .set { sample_ch }

Channel.fromPath('*.fq').map { file -> tuple(file.baseName, file) }
```

**4. `view`：查看 Channel 資料（除錯用）**

```groovy
my_channel.view()
```

* 可用於印出 Channel 中的資料內容。
* 也可傳入自定格式：

```groovy
my_channel.view { "Sample: $it" }
```

**5. `flatten`：攤平巢狀資料**

```groovy
Channel.of(['a', 'b'], ['c', 'd']).flatten().view()
```

* 將多層 list 展平成單一資料流：

  * 輸入：`[['a','b'],['c','d']]`
  * 輸出：`'a','b','c','d'`

---

### 🔹 與 Process 的互動

**輸入：**

```groovy
process align {
  input:
  tuple val(id), path(reads)

  script:
  """
  echo Processing $id with $reads
  """
}
```

**輸出：**

```groovy
output:
  path "*.bam" into bam_files
```

---

### 🔹 合併與展開

**合併多個 Channel：**

```groovy
Channel.merge(ch1, ch2)
```

**展開 tuple（`.collect()`）：**

```groovy
ch_input.map { id, file -> "${id}_new" }
```

---

`Channel` 是 Nextflow 的資料流基礎，可靈活操作以支援平行處理、條件分支與動態輸入。
若搭配 `emit:` 與 `into:` 可進一步模組化流程與改善可讀性。
