# MySQL JOIN 合併查詢教學

`JOIN` 用來從多個資料表中同時查詢資料。常見的 JOIN 類型包括 `INNER JOIN`、`LEFT JOIN`、`RIGHT JOIN` 與 `FULL JOIN`（MySQL 不直接支援 FULL JOIN，可用 UNION 模擬）。

## 資料表準備

```sql
-- 使用者資料表
CREATE TABLE users (
    id INT PRIMARY KEY,
    name VARCHAR(255)
);

-- 訂單資料表
CREATE TABLE orders (
    id INT PRIMARY KEY,
    user_id INT,
    product VARCHAR(255)
);

-- 插入範例資料
INSERT INTO users VALUES (1, 'Alice'), (2, 'Bob'), (3, 'Charlie');
INSERT INTO orders VALUES (101, 1, 'Book'), (102, 2, 'Pen'), (103, 2, 'Notebook');
```

## INNER JOIN

只回傳兩表皆有對應資料的列。

```sql
SELECT users.name, orders.product
FROM users
INNER JOIN orders ON users.id = orders.user_id;
```

## LEFT JOIN

保留左表所有資料，右表無對應則為 NULL。

```sql
SELECT users.name, orders.product
FROM users
LEFT JOIN orders ON users.id = orders.user_id;
```

## RIGHT JOIN

保留右表所有資料，左表無對應則為 NULL。

```sql
SELECT users.name, orders.product
FROM users
RIGHT JOIN orders ON users.id = orders.user_id;
```

## 模擬 FULL JOIN（使用 UNION）

```sql
SELECT users.name, orders.product
FROM users
LEFT JOIN orders ON users.id = orders.user_id
UNION
SELECT users.name, orders.product
FROM users
RIGHT JOIN orders ON users.id = orders.user_id;
```

---

`JOIN` 是進行關聯式資料查詢的重要工具，建議搭配 `EXPLAIN` 分析查詢效能，並確保資料表間具備適當的索引。
