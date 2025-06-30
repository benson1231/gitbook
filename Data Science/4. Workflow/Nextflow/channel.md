## Nextflow `Channel` 使用說明

在 Nextflow 中，`Channel` 是流程資料傳遞的核心。不同的 `process` 間透過 Channel 傳遞資料，以實現資料驅動的執行邏輯。

---

### 🔹 建立 Channel

#### 1. `fromPath`

從路徑建立檔案輸入通道。

```groovy
Channel.fromPath('data/*.fastq')
```

* **支援 glob 樣式**：如 `*.fq.gz`
* **傳出型別**：`path` 資料

#### 2. `from`

從 Iterable（如 List）建立通道。

```groovy
Channel.from(['sample1', 'sample2'])
```

* **多元素**：每個元素單獨 emit
* 可用於 `List`、`Set`、`Map` 等

#### 3. `of`

與 `from` 相同，但語法更簡潔。

```groovy
Channel.of('a', 'b', 'c')
```

* **立即指定多值**：可當作靜態輸入

#### 4. `value`

建立只包含一個資料項的通道。

```groovy
Channel.value('only_one')
```

* 與 `of('only_one')` 類似，但語義明確表示僅一值

#### 5. `empty`

建立空 Channel，預留用於後續 `.emit()` 動態填入。

```groovy
Channel.empty()
```

#### 6. `fromFilePairs`

建立成對檔案的資料流，常用於 R1/R2 成對資料。

```groovy
Channel.fromFilePairs('data/*_{R1,R2}.fastq.gz')
```

* 輸出為 `tuple(id, [file1, file2])`
* 常搭配 `input: tuple val(id), path(reads)`

#### 7. `fromPath(..., type: 'file')`

指定輸出為檔案（`file`）、目錄（`dir`）或混合。

---

### 🔹 Channel 操作方法

#### `map`

對通道內的每筆資料做轉換處理。

```groovy
Channel.from(['a','b']).map { it.toUpperCase() }
```

* 轉換值、加前綴/後綴、解析欄位等

#### `filter`

根據條件篩選通道資料。

```groovy
ch.filter { it.endsWith('.fastq') }
```

#### `flatten`

將嵌套 List 攤平為一層元素。

```groovy
Channel.of(['a', 'b'], ['c', 'd']).flatten()
```

* 輸入：`[['a','b'],['c','d']]`
* 輸出：`'a','b','c','d'`

#### `collect`

展開 tuple 結構供處理使用。

```groovy
ch.map { id, file -> "$id-processed" }
```

#### `set`

命名一個通道變數（DSL2 常用）。

```groovy
Channel.of('s1','s2').set { sample_ch }
```

#### `view`

除錯時顯示 channel 資料。

```groovy
my_ch.view()
my_ch.view { "Sample: $it" }
```

#### `merge`

合併多個 Channel 成一個。

```groovy
Channel.merge(ch1, ch2)
```

* **注意順序與資料型別一致性**

#### `mix`

交錯合併多個 Channel 的元素（較少用）。

```groovy
Channel.mix(ch1, ch2)
```

#### `combine`

兩個 channel 進行笛卡兒積組合。

```groovy
ch1.combine(ch2)
```

* 例如：`ch1 = [A,B], ch2 = [1,2]` → 輸出 `(A,1), (A,2), (B,1), (B,2)`

#### `ifEmpty`

為空時指定替代資料。

```groovy
ch.ifEmpty { Channel.of('default') }
```

#### `branch`

依據條件將資料分流。

```groovy
ch.branch {
  left: { it.startsWith('A') },
  right: { it.startsWith('B') }
}
```

* 回傳 `Map<String, Channel>`

#### `splitCsv`

將 CSV 逐列讀入並分欄成 tuple。

```groovy
Channel.fromPath('meta.csv').splitCsv(header:true)
```

* 搭配 `map` 可轉為 `(sample_id, file_path)`

#### `groupTuple`

將 tuple 通道依 key 分組。

```groovy
ch.groupTuple()
```

* 輸入：`(id, value)` → 輸出：`(id, [value1, value2, ...])`

#### `toList`

將所有 channel 內容收集成一個 list（結束時 emit）。

```groovy
ch.toList().view()
```

* 常見於 summary 或統計用途

---

### 🔹 與 Process 的互動

#### ➤ 作為輸入

```groovy
process align {
  input:
    tuple val(sample), path(reads)

  script:
    """
    echo Aligning $sample with $reads
    """
}
```

#### ➤ 作為輸出

```groovy
output:
  path '*.bam' into bam_ch
```

* 可配合 `into:` 將輸出導入下游 channel

---

### 🔹 實用技巧與補充

* `first()`：取得通道的第一筆資料（常用於 summary）
* `subscribe {}`：在腳本中以 callback 形式處理 channel 資料（偏 Groovy）
* `Channel.empty()` 搭配 `.emit:` 可用於流程動態建立輸出
* `Channel.fromPath().splitCsv()`：適用元資料處理
* `channel.toList()`：彙總輸出後統一處理（如 all BAM file）
* `channel.contains()`：在 workflow 中檢查是否包含某值

---

### 📌 小結

| 方法                | 說明                 |
| ----------------- | ------------------ |
| `fromPath()`      | 從檔案路徑建立通道          |
| `from()` / `of()` | 從 List 等集合建立資料流    |
| `value()`         | 建立只包含一個項目的通道       |
| `empty()`         | 建立空通道              |
| `map()`           | 資料轉換               |
| `filter()`        | 條件過濾               |
| `flatten()`       | 攤平成單層資料            |
| `view()`          | 印出資料流內容供除錯         |
| `merge()`         | 多通道合併為一            |
| `combine()`       | 兩通道資料配對組合（笛卡兒積）    |
| `ifEmpty()`       | 空通道給預設值            |
| `branch()`        | 分流到多通道             |
| `set()`           | 命名通道變數（DSL2）       |
| `splitCsv()`      | 讀取並解析 CSV 為欄位      |
| `groupTuple()`    | 將 tuple 通道依 key 分組 |
| `toList()`        | 將資料收集成 list 統一輸出   |

Nextflow 的 Channel 設計高度彈性，是串接、平行化流程的基石，學會各種建立與操作方法能大幅提高 pipeline 的模組化與擴充性。
