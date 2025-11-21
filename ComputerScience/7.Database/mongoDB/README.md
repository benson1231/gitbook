# MongoDB 介紹

[MongoDB](https://www.mongodb.com/) 是一種 NoSQL（非關聯式）資料庫，採用文件導向（document-oriented）設計，適合儲存結構彈性、變動頻繁的大量資料，廣泛用於現代 Web 應用與資料驅動的系統。

---

## ✅ MongoDB Compass 安裝

MongoDB 官方提供的 GUI 工具 —— **Compass**，能讓你視覺化瀏覽與操作資料庫。

### 安裝方式

👉 下載網址：[https://www.mongodb.com/try/download/compass](https://www.mongodb.com/try/download/compass)

1. 選擇作業系統（macOS / Windows / Linux）
2. 點選下載並安裝
3. 開啟 Compass 並填入連線字串，例如：

   ```
   mongodb://127.0.0.1:27017
   ```
4. 按下「Connect」，即可瀏覽本地 MongoDB

---

## 💻 MongoDB CLI 工具：mongosh

除了 GUI，也可使用 CLI 工具 `mongosh`（MongoDB Shell）與資料庫互動。

### 安裝（macOS）

```bash
brew install mongosh
```

### 基本用法

```bash
mongosh
```

連線成功後會看到提示字元：

```
test>
```

接著可執行：

```js
use mydb

db.users.insertOne({ name: "Alice", age: 30 })
db.users.find()
db.users.updateOne({ name: "Alice" }, { $set: { age: 31 } })
db.users.deleteOne({ name: "Alice" })
```

### 常用指令總覽

| 操作      | 指令範例                             |
| ------- | -------------------------------- |
| 切換資料庫   | `use mydb`                       |
| 新增文件    | `db.collection.insertOne({...})` |
| 查詢資料    | `db.collection.find()`           |
| 更新資料    | `db.collection.updateOne(...)`   |
| 刪除資料    | `db.collection.deleteOne(...)`   |
| 顯示所有資料庫 | `show dbs`                       |
| 顯示集合    | `show collections`               |

---

## 🧠 與傳統 SQL 資料庫比較

| 特性     | MongoDB (NoSQL)          | MySQL / PostgreSQL (SQL) |
| ------ | ------------------------ | ------------------------ |
| 資料結構   | JSON-like 文件（BSON 格式）    | 表格（Table）                |
| 資料單位   | 文件（Document）             | 列（Row）                   |
| 資料儲存方式 | 動態結構、可嵌套欄位               | 靜態欄位結構、正規化               |
| 語言     | MongoDB 查詢語法（類似 JSON）    | SQL 語法                   |
| 關聯處理   | 較弱，透過手動 join 或 `$lookup` | 支援 JOIN 等複雜關聯操作          |
| 可擴展性   | 水平擴展佳（sharding）          | 傳統垂直擴展為主                 |
| 適合應用   | 即時資料、物聯網、內容管理、原型設計等      | 金融、ERP、結構化資料密集型系統等       |

---

## 📦 MongoDB 資料結構概念

```
Database > Collection > Document
```

* **Database**：相當於 SQL 的資料庫
* **Collection**：相當於資料表（table）
* **Document**：一筆 JSON 資料，相當於一列（row）

### 文件範例

```json
{
  "name": "Alice",
  "age": 30,
  "email": "alice@example.com",
  "address": {
    "city": "Taipei",
    "zip": "100"
  }
}
```

---

## ✅ MongoDB 適合學習的對象與場景

* 前端/全端開發者（搭配 React、Vue、Node.js）
* 原型設計與快速部署專案
* 結構不固定、欄位動態變化的資料集
* 不需要大量 JOIN 或事務的系統

---

MongoDB 提供簡單、靈活的資料模型，是現代 Web 開發的常見選擇。搭配 Compass GUI 與 mongosh CLI，可快速上手並進行 CRUD 操作。若你習慣傳統 SQL，理解 MongoDB 的 document-based 架構將是邁向全端與雲端資料處理的重要一步。
