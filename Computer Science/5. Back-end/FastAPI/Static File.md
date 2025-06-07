# Static Files

在開發網站或應用程式時，常需要提供靜態資源，如圖片、CSS、JavaScript 檔案。FastAPI 提供 `StaticFiles` 模組來輕鬆託管這些靜態資源。

---

## 基本用法

靜態目錄：

```python
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

app = FastAPI()

# 將 ./static 資料夾掛載到 /static 路徑上
app.mount("/static", StaticFiles(directory="static"), name="static")
```

如果有以下檔案結構：

```
/static
  ├── style.css
  └── logo.png
```

即可透過以下網址存取：

* `http://localhost:8000/static/style.css`
* `http://localhost:8000/static/logo.png`

---

## 搭配 HTML 模板使用

靜態資源常與 Jinja2 模板結合，用於網頁渲染：

```html
<link rel="stylesheet" href="/static/style.css">
<img src="/static/logo.png">
```

---

## 實作示範：HTML + CSS

```python
from fastapi.responses import HTMLResponse

@app.get("/home", response_class=HTMLResponse)
def home():
    return """
    <html>
        <head>
            <link rel='stylesheet' href='/static/style.css'>
        </head>
        <body>
            <h1>Hello from static files!</h1>
            <img src='/static/logo.png'>
        </body>
    </html>
    """
```

---

## 小結

| 功能         | 用法                                   |
| ---------- | ------------------------------------ |
| 提供靜態資源     | `app.mount("/路徑", StaticFiles(...))` |
| 存取圖片、CSS 等 | `/static/檔名`                         |
| 結合 HTML 使用 | `<link>`、`<img>` 標籤引用                |

透過 `StaticFiles` 的整合，FastAPI 不僅能夠處理 API 邏輯，也能簡易支援網頁應用所需的靜態資源。
