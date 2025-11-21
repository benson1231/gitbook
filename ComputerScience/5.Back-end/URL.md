# Uniform Resource Locator (URL)

**last update: 2025-06-07**

`網址 (URL, Uniform Resource Locator)`是用來定位網路資源的標準格式。它包含多個部分，每個部分都有特定功能與意義。

---

## 範例網址

```
https://www.example.com:443/path/to/resource?search=keyword&page=2#section1
```

---

## 組成說明

| 部分                       | 說明                                        |
| ------------------------ | ----------------------------------------- |
| `https://`               | **協定（Protocol）**，決定通訊方式，常見為 HTTP 或 HTTPS  |
| `www.example.com`        | **主機名稱（Host）**，可為網域名稱或 IP 位址              |
| `:443`                   | **埠號（Port）**，可選，HTTPS 通常為 443，HTTP 為 80   |
| `/path/to/resource`      | **路徑（Path）**，指定伺服器上資源的位置                  |
| `?search=keyword&page=2` | **查詢參數（Query Parameters）**，鍵值對組合，用於傳遞額外資料 |
| `#section1`              | **片段識別（Fragment）**，錨點導向頁面內特定位置            |

---

## 查詢參數（Query Parameters）

查詢參數出現在 `?` 後面，格式為 `key=value`，多個參數以 `&` 分隔。

範例：

```
?q=python&page=3&lang=zh
```

等同於：

```json
{
  "q": "python",
  "page": "3",
  "lang": "zh"
}
```

---

## URL 實例解析

以 `http://127.0.0.1:8000/data?q=2&p=5` 為例：

```
http://127.0.0.1:8000/data?q=2&p=5
│     │         │       └───── 查詢參數（q=2, p=5）
│     │         └───────────── 路徑 /data
│     └──────────────────────── 埠號 8000
└────────────────────────────── 協定 http
```

---

網址的結構設計對於 Web 開發、API 規劃與 SEO 都具有重要意義。了解 URL 的組成，有助於正確解析、構建與測試各種網路應用。
