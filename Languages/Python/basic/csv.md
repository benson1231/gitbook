# 📚 Python CSV 全攻略（內建 `csv` 模組 + 常見實務）

> 這份教學聚焦 **Python 內建 `csv` 模組** 的正確用法與最佳實務，並補充少量 `pandas` 的進階技巧。涵蓋讀寫、方言（dialect）、引號/跳脫、編碼、巨量資料處理、常見陷阱與範例食譜。

---

## 🔰 基礎觀念

* **CSV**（Comma-Separated Values）是一種以分隔符號（預設逗號 `,`）分隔欄位的純文字格式。
* 各平台換行符可能不同（`\n`, `\r\n`）。
* 以 **Python 內建 `csv` 模組** 讀寫時，**開檔務必加上 `newline=''`**，以避免多出空白行（Windows 常見）。

```python
import csv
from pathlib import Path

p = Path('data.csv')
with p.open('w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['id', 'name', 'score'])
    writer.writerow([1, 'Alice', 95])
```

---

## 📖 讀取：`csv.reader` 與 `csv.DictReader`

### `csv.reader`

逐列回傳 **list**（每列一個 list）。

```python
import csv

with open('data.csv', 'r', newline='', encoding='utf-8') as f:
    reader = csv.reader(f)  # 預設 delimiter=','
    for row in reader:
        print(row)
```

常用參數：

* `delimiter=','`：分隔符（可改成 `\t` 處理 TSV）。
* `quotechar='"'`、`escapechar='\\'`：引號與跳脫。
* `skipinitialspace=True`：忽略分隔後的首個空白。

### `csv.DictReader`

每列回傳 **dict**（以首列或指定欄名當 key）。

```python
import csv

with open('data.csv', 'r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)  # 會讀取首列作為欄名
    for row in reader:
        print(row['name'], row.get('score'))
```

若檔案沒有標題列，可自訂欄名：

```python
with open('no_header.csv', 'r', newline='', encoding='utf-8') as f:
    fieldnames = ['id', 'name', 'score']
    reader = csv.DictReader(f, fieldnames=fieldnames)
    next(reader)  # 視情況略過第一列（若第一列是描述或無用資料）
```

### 讀取時的型別轉換

`csv` 讀進來都是字串，需自行轉型：

```python
with open('data.csv', 'r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    rows = []
    for r in reader:
        r['id'] = int(r['id'])
        r['score'] = float(r['score'])
        rows.append(r)
```

---

## ✍️ 寫入：`csv.writer` 與 `csv.DictWriter`

### `csv.writer`

```python
import csv

data = [
    ['id', 'name', 'note'],
    [1, 'Alice', 'He said: "Hi"'],  # 內含引號
    [2, 'Bob', 'Line1\nLine2']       # 內含換行
]

with open('out.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
    writer.writerows(data)
```

常用參數：

* `delimiter=','`：分隔符。
* `lineterminator='\n'`：列結尾（預設依平台）。
* `quoting=`：引號策略：

  * `csv.QUOTE_MINIMAL`（預設，必要時加引號）
  * `csv.QUOTE_ALL`（全欄位加引號）
  * `csv.QUOTE_NONNUMERIC`（非數值加引號，讀取時自動轉為 float）
  * `csv.QUOTE_NONE`（不加引號，通常需搭配 `escapechar`）
* `quotechar='"'`、`escapechar='\\'`

### `csv.DictWriter`

```python
import csv

rows = [
    {"id": 1, "name": "Alice", "score": 95},
    {"id": 2, "name": "Bob",   "score": 88},
]

with open('out_dict.csv', 'w', newline='', encoding='utf-8') as f:
    fieldnames = ['id', 'name', 'score']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
```

---

## 🗣️ Dialect（方言）與自動偵測

### 註冊自訂方言

```python
import csv

csv.register_dialect(
    'my_tsv',
    delimiter='\t',
    quoting=csv.QUOTE_MINIMAL,
    lineterminator='\n'
)

with open('data.tsv', 'r', newline='', encoding='utf-8') as f:
    reader = csv.reader(f, dialect='my_tsv')
```

### `csv.Sniffer` 自動偵測格式

```python
import csv

with open('unknown.csv', 'r', encoding='utf-8', newline='') as f:
    sample = f.read(2048)
    f.seek(0)
    dialect = csv.Sniffer().sniff(sample)
    has_header = csv.Sniffer().has_header(sample)

    reader = csv.reader(f, dialect)
    if has_header:
        headers = next(reader)
```

---

## 🌐 編碼（Encoding）與 Excel 相容

* 一般建議使用 **`utf-8`**。
* 若要給部分 Excel 版本（特別是 Windows/舊版）可靠辨識，寫檔時可用 **`utf-8-sig`** 以加入 BOM：

```python
with open('excel_friendly.csv', 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.writer(f)
    writer.writerow(['編號', '姓名'])
    writer.writerow([1, '王小明'])
```

* 若來源是 CP950/Big5 等本地編碼，讀取時需指定正確 `encoding`，並處理錯字：

```python
with open('legacy_big5.csv', 'r', newline='', encoding='cp950', errors='replace') as f:
    reader = csv.reader(f)
```

---

## 🚀 巨量檔案與效能技巧

* **逐行處理**：避免一次載入全部（特別是 `readlines()`），改用迭代器。
* **欄位過濾/映射**：在讀取迴圈中即時轉型與過濾，減少中間結構。
* **分塊寫出**：處理 N 筆就寫一次，避免巨大記憶體佔用。
* **壓縮檔**：可直接用 `gzip`/`bz2` 搭配檔案物件：

```python
import csv, gzip

with gzip.open('data.csv.gz', 'rt', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        ...
```

---

## 🧪 常見陷阱與排錯

1. **Windows 寫檔多出空行**：忘了 `newline=''`。
2. **欄位內含逗號/引號/換行**：需設定 `quoting`（預設即可處理），或搭配 `escapechar`。
3. **型別全是字串**：自行轉型（`int`, `float`, `datetime`）。
4. **雜訊標頭/註解列**：讀取時先 `next(reader)` 跳過或預處理。
5. **分隔符不是逗號**：指定 `delimiter` 或用 `Sniffer`。

---

## 🍱 實用食譜（Recipes）

### 1) 篩選列並輸出新檔

```python
import csv

with open('sales.csv', 'r', newline='', encoding='utf-8') as fin, \
     open('high_value.csv', 'w', newline='', encoding='utf-8') as fout:

    reader = csv.DictReader(fin)
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()

    for r in reader:
        if float(r['amount']) >= 1000:
            writer.writerow(r)
```

### 2) 轉置（rows ↔ columns）

```python
import csv

with open('in.csv', 'r', newline='', encoding='utf-8') as f:
    rows = list(csv.reader(f))

transposed = list(zip(*rows))

with open('out_transposed.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerows(transposed)
```

### 3) 追加新欄位（依既有欄計算）

```python
import csv

with open('scores.csv', 'r', newline='', encoding='utf-8') as fin, \
     open('scores_plus.csv', 'w', newline='', encoding='utf-8') as fout:

    reader = csv.DictReader(fin)
    fieldnames = reader.fieldnames + ['passed']
    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()

    for r in reader:
        r['passed'] = 'Y' if float(r['score']) >= 60 else 'N'
        writer.writerow(r)
```

### 4) 合併多個 CSV（同欄位）

```python
import csv, glob

files = glob.glob('parts/part_*.csv')
with open('merged.csv', 'w', newline='', encoding='utf-8') as fout:
    writer = None
    for path in files:
        with open(path, 'r', newline='', encoding='utf-8') as fin:
            reader = csv.DictReader(fin)
            if writer is None:
                writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
                writer.writeheader()
            for row in reader:
                writer.writerow(row)
```

### 5) 自訂分隔符（TSV / 管線符）

```python
import csv

with open('data.tsv', 'r', newline='', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        ...

with open('pipe.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='|')
    writer.writerow(['a', 'b', 'c'])
```

---

## 🧰 與 `pandas` 的快速對照

* **讀**：

```python
import pandas as pd

df = pd.read_csv('data.csv', encoding='utf-8')
# 常用參數：usecols, dtype, parse_dates, na_values, chunksize
```

* **寫**：

```python
df.to_csv('out.csv', index=False, encoding='utf-8-sig')
```

* **大量資料**：`chunksize` 分塊讀取、指定 `dtype` 降低記憶體。

---

## ✅ 重點總結

1. 讀寫 CSV 時 **一定用 `newline=''`**；編碼建議 `utf-8`/`utf-8-sig`（對 Excel 友善）。
2. 欄位內含逗號、引號、換行時，調整 `quoting`/`escapechar`；或保留預設讓模組自動處理。
3. `DictReader/DictWriter` 能以欄名取值，更可讀。
4. 巨量檔案：逐行處理、分塊寫出、必要時搭配壓縮與 `pandas`。
5. 善用 `Sniffer` 與 Dialect 對付未知或非標準 CSV。
