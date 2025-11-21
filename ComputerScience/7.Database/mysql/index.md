# MySQL 索引（INDEX）教學

索引（Index）是用來加快資料表查詢效率的資料結構，尤其在大量資料或複雜條件查詢中能顯著提升效能。

## 為什麼需要索引？

* 加速 `SELECT` 查詢速度。
* 提高 `JOIN`、`WHERE` 條件與 `ORDER BY` 操作效能。
* 強制欄位唯一（`UNIQUE` 索引）。

> 注意：雖然索引能提升查詢效能，但也會增加 `INSERT`、`UPDATE`、`DELETE` 操作的成本，並佔用額外的儲存空間。

## 常見索引類型

| 類型          | 說明                                       |
| ----------- | ---------------------------------------- |
| PRIMARY KEY | 主鍵索引，自動建立唯一且非空索引                         |
| UNIQUE      | 唯一索引，避免重複值                               |
| INDEX / KEY | 一般索引，可加速查詢但允許重複值                         |
| FULLTEXT    | 全文索引，用於文字搜尋（僅支援 `MyISAM` 或 `InnoDB` 新版本） |
| SPATIAL     | 空間索引，用於地理資訊欄位                            |

### 索引的資料結構：B-Tree

MySQL 中大部分的索引預設使用 B-Tree 結構。

* **B-Tree（Balanced Tree）**：適用於範圍查詢（例如 `BETWEEN`、`>`, `<`）與排序操作。
* **Hash Index**：主要由 Memory 引擎支援，僅適合等值查詢（不支援範圍查詢）。

## 建立索引方式

### 1. 使用 `CREATE INDEX`

```sql
CREATE INDEX idx_emp_name ON employees(emp_name);
```

### 2. 建立唯一索引

```sql
CREATE UNIQUE INDEX idx_email_unique ON users(email);
```

### 3. 建立複合索引（多欄位）

```sql
CREATE INDEX idx_name_age ON employees(emp_name, age);
```

## 使用 `ALTER TABLE` 新增索引

```sql
ALTER TABLE employees ADD INDEX idx_name (emp_name);
ALTER TABLE users ADD UNIQUE INDEX idx_email (email);
ALTER TABLE orders ADD INDEX idx_user_product (user_id, product_id);
```

## 在建立資料表時加入索引

```sql
CREATE TABLE users (
    id INT PRIMARY KEY,
    name VARCHAR(100),
    email VARCHAR(100),
    INDEX idx_name (name),
    UNIQUE KEY idx_email (email)
);
```

## 刪除索引

```sql
DROP INDEX idx_name ON users;
```

## 查詢資料表索引資訊

```sql
SHOW INDEX FROM users;
```

---

適當的索引設計可大幅提升資料庫查詢效能，建議搭配 `EXPLAIN` 分析查詢計畫，以找出最佳化的欄位索引組合。

> 📌 小技巧：避免對選擇性（selectivity）太低的欄位建立索引，例如性別（M/F）或布林值，否則可能反而降低效能。
