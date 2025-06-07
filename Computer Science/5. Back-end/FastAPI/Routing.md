# Routing

在 Web 開發中，\*\*路由（Routing）\*\*指的是將使用者輸入的網址路徑，對應到後端系統中的特定功能或資料處理邏輯。

---

## 路由的基本概念

路由通常定義一組 URL 模式與對應的處理函數（或稱處理器、controller）。例如：

```python
@app.get("/")
def home():
    return {"message": "這是首頁"}
```

---

## 靜態路由（Static Routing）

靜態路由是 URL 結構固定、不可變的路徑。

範例：

```python
@app.get("/about")
def about():
    return {"info": "這是關於我們的頁面"}
```

呼叫：`http://localhost:8000/about`

---

## 動態路由（Dynamic Routing）

動態路由允許路徑中包含變數部分，可對應不同參數內容。

### FastAPI 範例：

```python
@app.get("/users/{user_id}")
def get_user(user_id: int):
    return {"user_id": user_id}
```

呼叫：`http://localhost:8000/users/42` → 回傳 `{"user_id": 42}`

### 多個動態參數：

```python
@app.get("/books/{category}/{book_id}")
def book_detail(category: str, book_id: int):
    return {"category": category, "id": book_id}
```

---

## 路由參數的型別註記

FastAPI 支援明確的型別註記與驗證，例如：

```python
@app.get("/item/{id}")
def get_item(id: int):
    return {"result": id * 2}
```

如果使用者輸入非整數（如 `/item/abc`），將自動回傳錯誤訊息。

---

## 查詢參數（Query Parameters）與路由配合

```python
@app.get("/search")
def search(q: str, page: int = 1):
    return {"q": q, "page": page}
```

呼叫：`/search?q=fastapi&page=2`

---

## 小結

| 類型   | 範例 URL                    | 使用說明            |
| ---- | ------------------------- | --------------- |
| 靜態路由 | `/about`                  | 對應固定功能，如首頁、關於頁面 |
| 動態路由 | `/user/{user_id}`         | 使用者變數輸入、依參數動態處理 |
| 查詢參數 | `/search?q=python&page=2` | 用於傳遞篩選、查詢、分頁等條件 |

理解靜態與動態路由的使用方式，能幫助你設計清晰、有彈性的 API 系統。
