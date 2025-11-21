# requests

`requests` 是 Python 最常用的第三方套件之一，用於處理 HTTP 請求，包含 GET、POST、PUT、DELETE 等常見操作，並可與 JSON、API、認證等功能整合使用。

---

## 一、安裝 requests 套件

```bash
pip install requests
```

---

## 二、基本用法：GET 請求

```python
import requests

response = requests.get('https://jsonplaceholder.typicode.com/posts/1')

print(response.status_code)  # HTTP 狀態碼
print(response.headers['Content-Type'])  # 回應標頭
print(response.json())  # 回應 JSON 內容（自動轉換為 dict）
```

---

## 三、POST 請求（含 JSON 資料傳送）

```python
data = {'title': 'foo', 'body': 'bar', 'userId': 1}
response = requests.post(
    'https://jsonplaceholder.typicode.com/posts',
    json=data
)

print(response.status_code)
print(response.json())
```

---

## 四、常見 HTTP 方法

| 方法     | 說明    | 範例語法                            |
| ------ | ----- | ------------------------------- |
| GET    | 取得資源  | `requests.get(url)`             |
| POST   | 建立新資源 | `requests.post(url, data=...)`  |
| PUT    | 完整更新  | `requests.put(url, data=...)`   |
| PATCH  | 部分更新  | `requests.patch(url, data=...)` |
| DELETE | 刪除資源  | `requests.delete(url)`          |

---

## 五、加入 Headers、Params、Timeout

```python
headers = {'Authorization': 'Bearer my-token'}
params = {'q': 'python'}

response = requests.get(
    'https://example.com/api',
    headers=headers,
    params=params,
    timeout=5  # 最長等待秒數
)
```

---

## 六、處理錯誤與例外

```python
try:
    r = requests.get('https://example.com', timeout=3)
    r.raise_for_status()  # 檢查是否為 4xx 或 5xx 錯誤
except requests.exceptions.RequestException as e:
    print("請求發生錯誤：", e)
```

---

## 七、應用場景

* 串接第三方 API（如天氣、翻譯、GitHub）
* 上傳表單、檔案或圖片（搭配 `files` 參數）
* 與 Flask/Django 建構的 RESTful API 溝通
* 搭配 JSON 處理資料交換

---

`requests` 是進行網路通訊與資料擷取的重要工具，簡潔易用、功能完整，適合幾乎所有 Python 網路應用情境。
