# Python JSON 模組教學

## 一、基本介紹

`json` 是 Python 內建的模組，用於 **JSON (JavaScript Object Notation)** 的讀寫與轉換。
它常被用於：

* Web API 資料交換（如 RESTful API）
* 檔案儲存設定或結果
* Python 與 JavaScript 的資料交互

Python 與 JSON 的資料型別對應如下：

| Python 型別   | JSON 型別 |
| ----------- | ------- |
| dict        | object  |
| list, tuple | array   |
| str         | string  |
| int, float  | number  |
| True        | true    |
| False       | false   |
| None        | null    |

---

## 二、常用方法

### 1. `json.dumps()` – 將 Python 物件轉成 JSON 字串

```python
import json

data = {"name": "Alice", "age": 25, "city": "Taipei"}
json_str = json.dumps(data)
print(json_str)
# 輸出: '{"name": "Alice", "age": 25, "city": "Taipei"}'
```

#### 常見參數：

* `indent=4`：輸出格式化 JSON 字串。
* `sort_keys=True`：依鍵名排序。
* `ensure_ascii=False`：允許輸出中文。

```python
print(json.dumps(data, indent=4, ensure_ascii=False))
```

---

### 2. `json.loads()` – 將 JSON 字串轉成 Python 物件

```python
json_str = '{"name": "Bob", "age": 30}'
data = json.loads(json_str)
print(data["name"])  # Bob
```

---

### 3. `json.dump()` – 將 Python 物件寫入 JSON 檔

```python
with open('data.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=4, ensure_ascii=False)
```

> 適合直接輸出到檔案。

---

### 4. `json.load()` – 從檔案讀取 JSON

```python
with open('data.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
print(data)
```

> 常用於讀取設定檔、API 回應等。

---

## 三、錯誤處理

當 JSON 格式錯誤時，會拋出 `json.JSONDecodeError`。

```python
from json import loads, JSONDecodeError

try:
    data = loads('{name: "Bob"}')  # 缺少引號，錯誤
except JSONDecodeError as e:
    print("無效的 JSON 格式:", e)
```

---

## 四、實務應用範例

### 1️⃣ 儲存程式設定

```python
config = {"version": 1.2, "debug": True, "output": "result.txt"}
with open('config.json', 'w', encoding='utf-8') as f:
    json.dump(config, f, indent=2)
```

### 2️⃣ 從 API 解析資料

```python
import requests, json

response = requests.get('https://api.github.com/users/octocat')
data = json.loads(response.text)
print(data['login'])
```

### 3️⃣ 自訂編碼/解碼（進階）

若物件非原生 JSON 格式，可用 `default` 參數自訂轉換：

```python
import json, datetime

def encoder(obj):
    if isinstance(obj, datetime.datetime):
        return obj.isoformat()
    raise TypeError(f'Type {type(obj)} not serializable')

now = {"time": datetime.datetime.now()}
print(json.dumps(now, default=encoder, indent=2))
```

---

## 五、總結

| 函式             | 功能               | 常用場景     |
| -------------- | ---------------- | -------- |
| `json.dumps()` | Python → JSON 字串 | 傳輸、日誌    |
| `json.loads()` | JSON 字串 → Python | API 回應解析 |
| `json.dump()`  | Python → JSON 檔案 | 儲存設定或結果  |
| `json.load()`  | JSON 檔案 → Python | 讀取設定檔    |

JSON 模組是 Python 與外部世界溝通的核心之一，熟練掌握它能顯著提升資料處理與系統整合能力。
