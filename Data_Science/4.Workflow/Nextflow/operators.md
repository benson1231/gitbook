# Nextflow Operators 教學筆記

Nextflow 的 operators 是用來 **操作與轉換 channels** 的方法。它們支援 **Dataflow programming** 的風格，能將輸入資料流轉換為新的 channel。除了 `set` 與 `subscribe`，幾乎所有 operators 都會輸出新的 channel，因此能夠進行 **chaining**。

Operators 主要可分為七大類：

* **Filtering operators**
* **Transforming operators**
* **Splitting operators**
* **Combining operators**
* **Forking operators**
* **Maths operators**
* **Other operators**

---

## 1. 基本範例

最常見的 `map` operator：

```groovy
nums = Channel.of(1, 2, 3, 4)
square = nums.map { it -> it * it }
square.view()
```

輸出：

```
1
4
9
16
```

operators 可以鏈接：

```groovy
Channel
    .of(1, 2, 3, 4)
    .map { it * it }
    .view()
```

---

## 2. 常用 Operators

### 2.1 `view()`

將 channel 內容輸出到 console。

```groovy
Channel.of('foo', 'bar').view()
// 輸出:
foo
bar
```

也可以自訂格式：

```groovy
Channel.of('foo', 'bar').view { "- $it" }
```

---

### 2.2 `map()`

轉換每個元素，回傳新的 channel。

```groovy
Channel.of('hello', 'world')
    .map { it.reverse() }
    .view()
```

輸出：

```
olleh
dlrow
```

也可輸出 tuple：

```groovy
Channel.of('hello', 'world')
    .map { w -> [w, w.size()] }
    .view()
```

輸出：

```
[hello, 5]
[world, 5]
```

---

### 2.3 `mix()`

合併多個 channels：

```groovy
c1 = Channel.of(1, 2)
c2 = Channel.of('a', 'b')

c1.mix(c2).view()
```

輸出 (順序不保證)：

```
1
2
a
b
```

---

### 2.4 `flatten()`

展開 tuple 或 list：

```groovy
Channel.of([1, 2], [3, 4])
    .flatten()
    .view()
```

輸出：

```
1
2
3
4
```

---

### 2.5 `collect()`

收集所有輸出為一個 list（**value channel**）：

```groovy
Channel.of(1, 2, 3)
    .collect()
    .view()
```

輸出：

```
[1, 2, 3]
```

---

### 2.6 `groupTuple()`

依 key 分組：

```groovy
Channel.of([1,'A'], [1,'B'], [2,'C'])
    .groupTuple()
    .view()
```

輸出：

```
[1, [A, B]]
[2, [C]]
```

---

### 2.7 `join()`

以 key (tuple 第一個元素) 合併：

```groovy
left  = Channel.of(['X',1], ['Y',2])
right = Channel.of(['Y',5], ['X',4])

left.join(right).view()
```

輸出：

```
[Y, 2, 5]
[X, 1, 4]
```

---

### 2.8 `branch()`

根據條件分流：

```groovy
Channel.of(1, 2, 30)
    .branch {
        small: it < 10
        large: it >= 10
    }
    .set { result }

result.small.view { "$it is small" }
result.large.view { "$it is large" }
```

輸出：

```
1 is small
2 is small
30 is large
```

---

## 3. 文字檔案處理

### 3.1 `splitText()`

將文字檔依行數切分。

```groovy
Channel.fromPath('data/meta/random.txt')
    .splitText(by: 2)
    .view()
```

也可轉換：

```groovy
Channel.fromPath('data/meta/random.txt')
    .splitText(by: 2) { it.toUpperCase() }
    .view()
```

---

### 3.2 `splitCsv()`

解析 CSV：

```groovy
Channel.fromPath('data/meta/patients_1.csv')
    .splitCsv(header: true)
    .view { row -> "${row.patient_id}, ${row.num_samples}" }
```

可處理多個 CSV：

```groovy
Channel.fromPath('data/meta/patients_*.csv')
    .splitCsv(header: true)
    .view()
```

---

### 3.3 `splitJson()`

解析 JSON：

```groovy
Channel.of('[{"name":"Bob"},{"name":"Alice"}]')
    .splitJson()
    .view()
```

輸出：

```
[name:Bob]
[name:Alice]
```

也可直接從檔案：

```groovy
Channel.fromPath('file.json')
    .splitJson()
    .view()
```

---

## 4. 小結

* `view()`：直接輸出
* `map()`：逐項轉換
* `mix()`：合併多個 channels
* `flatten()`：展開 list/tuple
* `collect()`：一次性收集
* `groupTuple()`：依 key 分組
* `join()`：兩 channel 依 key 合併
* `branch()`：依條件分流
* `splitText()` / `splitCsv()` / `splitJson()`：處理文字檔

---

📌 建議：Operators 常與 **channel factories** 搭配使用（例如 `fromPath` → `.splitCsv()` → `.map()`），能夠快速建構彈性強的 workflow。
