# JSON in Python

JSON（JavaScript Object Notation）是一種輕量級的資料交換格式，廣泛應用於 API 資料傳輸與設定檔儲存。Python 提供內建的 `json` 模組來讀寫與轉換 JSON 資料。

---

## 一、基本語法與資料型態對應

| JSON 類型      | Python 類型    |
| ------------ | ------------ |
| object       | dict         |
| array        | list         |
| string       | str          |
| number       | int / float  |
| true / false | True / False |
| null         | None         |

---

## 二、讀取 JSON 檔案（json.load）

```python
import json

with open('data.json', 'r', encoding='utf-8') as f:
    data = json.load(f)
    print(data)
```

---

## 三、寫入 JSON 檔案（json.dump）

```python
person = {'name': 'Alice', 'age': 30, 'is_student': False}

with open('output.json', 'w', encoding='utf-8') as f:
    json.dump(person, f, ensure_ascii=False, indent=2)
```

### `indent` 參數可讓輸出更具可讀性

### `ensure_ascii=False` 支援中文輸出

---

## 四、JSON 字串處理（load vs loads）

```python
# 將 JSON 字串轉成 Python 物件
json_str = '{"x": 10, "y": [1, 2, 3]}'
data = json.loads(json_str)

# 將 Python 物件轉為 JSON 字串
json_text = json.dumps(data, indent=2)
print(json_text)
```

---

## 五、錯誤處理與常見問題

| 問題                 | 原因與解法                         |
| ------------------ | ----------------------------- |
| UnicodeEncodeError | 輸出中文時未設定 `ensure_ascii=False` |
| TypeError          | 資料中有無法序列化的物件（如 datetime、物件）   |
| JSONDecodeError    | 輸入 JSON 格式不正確                 |

---

## 六、應用情境

* 儲存程式設定檔、使用者偏好（config.json）
* 資料交換（API 請求與回應）
* 輕量級資料儲存格式（如取代 CSV 的結構性需求）

---

熟悉 JSON 的讀寫與轉換是現代程式開發的重要技能，善用 `json` 模組能讓你輕鬆處理結構化資料與串接各類服務。
