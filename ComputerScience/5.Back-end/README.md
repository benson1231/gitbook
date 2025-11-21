# Back-end

後端是指網站或應用程式中負責處理資料、邏輯、使用者認證、伺服器設定、API 提供等運算的部分，通常與資料庫整合並回應前端的請求。

---

## 核心概念

* **伺服器（Server）**：處理來自用戶端（client）的請求並提供對應回應。
* **路由（Routing）**：定義 URL 對應的功能處理程式。
* **資料庫（Database）**：儲存應用程式所需的資料。
* **API（Application Programming Interface）**：提供前端與後端之間的資料交換介面，常使用 REST 或 GraphQL 架構。

---

## 常見後端語言與框架

| 語言         | 框架 / 工具                  | 特點               |
| ---------- | ------------------------ | ---------------- |
| Python     | Flask / FastAPI / Django | 簡潔易讀，適合資料處理與原型開發 |
| JavaScript | Node.js + Express        | 單一語言前後端整合        |
| Java       | Spring Boot              | 嚴謹、適合大型企業系統      |
| PHP        | Laravel / Symfony        | 經典網頁後端語言，資源豐富    |
| Go         | Gin / Echo               | 高效能，適合微服務與高併發    |

---

## FastAPI 範例（Python）

```python
from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def read_root():
    return {"message": "Hello, backend!"}
```

---

## RESTful API 介紹

REST（Representational State Transfer）是 Web API 的主流架構，其設計風格如下：

| HTTP 方法 | 用途   | 範例 URI     |
| ------- | ---- | ---------- |
| GET     | 取得資料 | `/items/1` |
| POST    | 新增資料 | `/items/`  |
| PUT     | 更新資料 | `/items/1` |
| DELETE  | 刪除資料 | `/items/1` |

---

## 資料格式

後端與前端溝通時常使用 JSON 格式：

```json
{
  "id": 1,
  "name": "sample",
  "status": "active"
}
```

---

## 常見後端任務

* 使用者登入與驗證（Authentication / Authorization）
* 表單處理與檢查
* 資料查詢與 CRUD 操作
* 與外部 API 整合（如金流、地圖、AI 模型）
* 日誌紀錄與錯誤處理

---

後端開發為一個網站或服務的核心邏輯提供支撐，需與前端緊密合作，並確保系統穩定與資料安全。選擇合適的語言與框架能讓開發事半功倍。
