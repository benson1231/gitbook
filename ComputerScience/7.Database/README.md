# Database

資料庫是資訊系統的核心之一，負責有系統地儲存、查詢與管理資料。在電腦科學中，資料庫廣泛應用於網站、應用程式、雲端服務、人工智慧、大數據等領域。

---

## 一、資料庫分類

### 1. 關聯式資料庫（Relational Database, RDBMS）

* 基於表格（Table）與主鍵（Primary Key）設計
* 使用 SQL（Structured Query Language）查詢與操作
* 適用場景：事務性系統、資料一致性要求高的應用
* 代表系統：

  * MySQL / MariaDB
  * PostgreSQL
  * SQLite
  * Oracle DB
  * Microsoft SQL Server

### 2. 非關聯式資料庫（NoSQL）

* 不以表格為核心，適合非結構化或半結構化資料
* 類型：

  * 文件型（Document-Based）：MongoDB、CouchDB
  * 鍵值型（Key-Value）：Redis、DynamoDB
  * 廣欄型（Wide Column）：Cassandra、HBase
  * 圖形型（Graph-Based）：Neo4j、ArangoDB
* 優點：可擴展性高、靈活性強、效能佳

---

## 二、關鍵概念

| 概念                   | 說明                             |
| -------------------- | ------------------------------ |
| 資料模型                 | 描述資料如何結構與關聯（例：ER 模型、階層模型、文件模型） |
| 正規化（Normalization）   | 拆解冗餘資料，提高一致性與查詢效率              |
| 交易（Transaction）      | 一組資料操作的單位，需滿足 ACID 特性          |
| 索引（Index）            | 加速查詢速度，類似書籍的目錄結構               |
| 查詢語言（Query Language） | 最常見為 SQL，也有 NoSQL 系統的專屬查詢語法    |

---

## 三、SQL 語法簡介

```sql
-- 建立資料表
CREATE TABLE users (
  id INT PRIMARY KEY,
  name TEXT,
  email TEXT
);

-- 查詢
SELECT * FROM users WHERE name LIKE 'A%';

-- 更新
UPDATE users SET email='new@example.com' WHERE id=1;

-- 刪除
DELETE FROM users WHERE id=1;
```

---

## 四、現代資料庫應用

* **雲端資料庫**：Amazon RDS、Google Cloud SQL、Firebase
* **時序資料庫**：InfluxDB、TimescaleDB（用於 IoT 或監控）
* **全文搜尋資料庫**：Elasticsearch（支援模糊查詢與分詞）
* **混合式資料庫（HTAP）**：同時支援 OLTP 與 OLAP，如 TiDB

---

## 五、資料庫選擇建議

| 情境                | 建議資料庫               |
| ----------------- | ------------------- |
| 高度結構化資料 + ACID 保證 | PostgreSQL, MySQL   |
| 文件型 JSON 資料       | MongoDB             |
| 快速鍵值存取 / 快取       | Redis               |
| 巨量串流 / 可擴展性要求高    | Cassandra, DynamoDB |
| 關係圖譜 / 網路結構       | Neo4j               |

---

資料庫技術是每個電腦科學領域專業人士必備的核心知識，無論是資料科學、軟體工程還是系統設計，其設計原則與選型考量都直接影響系統的效能與可靠性。
