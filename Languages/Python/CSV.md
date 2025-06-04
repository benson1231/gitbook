# CSV in Python

CSV（Comma-Separated Values）是常見的純文字資料格式，廣泛應用於資料儲存與交換。Python 提供數種方式操作 CSV 檔案，最常用的是標準函式庫 `csv` 模組與 `pandas` 套件。

---

## 一、使用 `csv` 模組操作 CSV

### 讀取 CSV：

```python
import csv

with open('data.csv', newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        print(row)
```

### 寫入 CSV：

```python
with open('output.csv', 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Name', 'Age'])
    writer.writerow(['Alice', 25])
```

### 使用 DictReader / DictWriter：

```python
# 讀取成字典
with open('data.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        print(row['Name'], row['Age'])

# 寫入成字典
with open('output.csv', 'w', newline='') as csvfile:
    fieldnames = ['Name', 'Age']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerow({'Name': 'Bob', 'Age': 30})
```

---

## 二、使用 Pandas 操作 CSV

### 讀取 CSV：

```python
import pandas as pd

df = pd.read_csv('data.csv')
print(df.head())
```

### 寫入 CSV：

```python
df.to_csv('output.csv', index=False)
```

### 常見操作：

```python
# 篩選欄位
print(df[['Name', 'Age']])

# 加總數值欄位
print(df['Salary'].sum())

# 排序
print(df.sort_values(by='Age'))
```

---

## 三、常見注意事項

| 問題類型   | 解決方法                            |
| ------ | ------------------------------- |
| 編碼錯誤   | 指定 `encoding='utf-8'` 或 `big5`  |
| 換行問題   | 使用 `newline=''` 避免空行            |
| 欄位名稱重複 | `DictReader` 可解決對應問題            |
| 欄位順序   | `DictWriter` 寫入時指定 `fieldnames` |

---

## 四、適用情境比較

| 工具       | 優點           | 適合情境             |
| -------- | ------------ | ---------------- |
| `csv`    | 內建模組、無須額外安裝  | 小型 CSV、簡易格式轉換    |
| `pandas` | 功能強大、資料處理效率高 | 資料分析、大型資料集、結構化運算 |

---

CSV 操作是資料處理流程中的關鍵步驟，無論是初步匯入資料或將處理結果輸出，善用 `csv` 與 `pandas` 模組能提升開發效率與可維護性。
