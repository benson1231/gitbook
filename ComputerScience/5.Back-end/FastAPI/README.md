# FastAPI

**last update: 2025-06-07**

FastAPI 是一個現代化、快速（高效能）的 Python Web 框架，專為建立 RESTful API 而設計，支援自動生成文件（Swagger UI），並能進行型別驗證整合。

---

## 安裝 FastAPI 與伺服器

使用 pip 安裝 `FastAPI` 及 ASGI 伺服器 `uvicorn`：

```bash
pip install fastapi
pip install "uvicorn[standard]"
```

---

## 建立主程式 `main.py`

以下為簡單的 FastAPI 範例：

```python
from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def index():
    return {"message": "Hello World", "year": 2025}
```

---

## 啟動伺服器

使用 `uvicorn` 啟動 FastAPI 應用：

```bash
uvicorn main:app --reload

# 自訂埠號
uvicorn main:app --reload --port 3000
```

* `main:app`：`main.py` 檔案中的 `app` 實例
* `--reload`：啟用自動重新加載功能（適合開發階段）

---

## 使用說明

啟動後，可在以下網址瀏覽：

* API 主頁（Hello World）：`http://127.0.0.1:8000/`
* 傳入參數（例如 `/data?q=2&p=5`）：`http://127.0.0.1:8000/data?q=2&p=5`
* Swagger 自動文件介面：`http://127.0.0.1:8000/docs`

---

## 結語

`FastAPI` 是目前最受歡迎的 Python API 框架之一，具備簡潔語法、強型別驗證與自動文件支援，適合快速開發現代 Web API。
