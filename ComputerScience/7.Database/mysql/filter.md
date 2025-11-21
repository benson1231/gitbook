# MySQL 條件查詢指令教學（WHERE）

`WHERE` 子句用來根據條件篩選資料。以下範例說明常見的條件運算與應用方式。

## 資料表準備

```sql
CREATE TABLE mytable(
    id INT PRIMARY KEY AUTO_INCREMENT,
    age INT NOT NULL DEFAULT 18,
    name VARCHAR(255) NOT NULL
);

INSERT INTO mytable(age, name) VALUES (18, 'benson');
INSERT INTO mytable(age, name) VALUES (20, 'evan');
INSERT INTO mytable(age, name) VALUES (25, 'alice');
INSERT INTO mytable(age, name) VALUES (30, 'bob');
```

## 使用 `WHERE` 進行條件查詢

### 1. 單一條件

```sql
-- 查詢 age = 18 的資料
SELECT * FROM mytable WHERE age = 18;
```

### 2. 不等於條件（<>）

```sql
-- 查詢 age 不等於 18 的資料
SELECT * FROM mytable WHERE age <> 18;
```

### 3. 其他常見比較運算子

```sql
-- age 大於 20
SELECT * FROM mytable WHERE age > 20;

-- age 小於或等於 25
SELECT * FROM mytable WHERE age <= 25;
```

### 4. 多重條件

```sql
-- age 大於 20 且 name 為 'bob'
SELECT * FROM mytable WHERE age > 20 AND name = 'bob';

-- age 小於 30 或 name 為 'alice'
SELECT * FROM mytable WHERE age < 30 OR name = 'alice';
```

### 5. 模糊比對 (`LIKE`)

```sql
-- name 開頭為 'b'
SELECT * FROM mytable WHERE name LIKE 'b%';

-- name 結尾為 'n'
SELECT * FROM mytable WHERE name LIKE '%n';
```

---

`WHERE` 子句是資料篩選的基礎，搭配不同運算子可以實現各種查詢條件。
