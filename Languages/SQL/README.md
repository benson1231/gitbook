# SQL 指令教學

#### [SQL指令](https://www.codecademy.com/article/sql-commands)

#### [下載PostgreSQL](https://www.codecademy.com/article/installing-and-using-postgresql-locally)

#### [sqliteviz(線上操作)](https://sqliteviz.com/app/#/)

#### [papaya教室公開範例](https://papayaclassroom.notion.site/SQL-eb3a9ce2c9404f518674d885c5a789a5)

---

## 🧭 簡介

SQL（Structured Query Language）是用於與關聯式資料庫互動的語言。透過 SQL，你可以查詢（query）、插入（insert）、更新（update）與刪除（delete）資料，也可以建立與修改資料表結構。

SQL 指令可分為四大類：

* **DDL (Data Definition Language)**：定義資料庫結構，例如 `CREATE`、`ALTER`、`DROP`。
* **DML (Data Manipulation Language)**：操作資料內容，例如 `SELECT`、`INSERT`、`UPDATE`、`DELETE`。
* **DCL (Data Control Language)**：管理使用者權限，例如 `GRANT`、`REVOKE`。
* **TCL (Transaction Control Language)**：控制交易，例如 `COMMIT`、`ROLLBACK`。

---

## 📘 常見 SQL 指令

### 1️⃣ SELECT

從資料表中擷取資料。

```sql
SELECT column1, column2 FROM table_name;
```

可使用 `*` 選取所有欄位：

```sql
SELECT * FROM employees;
```

加入條件：

```sql
SELECT name, age FROM employees WHERE age > 30;
```

---

### 2️⃣ DISTINCT — 去除重複資料

`DISTINCT` 用於刪除查詢結果中重複的值。

```sql
SELECT DISTINCT department FROM employees;
```

如果想查詢部門與職稱的唯一組合：

```sql
SELECT DISTINCT department, job_title FROM employees;
```

---

### 3️⃣ INSERT

向資料表新增資料列。

```sql
INSERT INTO employees (name, age, department) VALUES ('Alice', 29, 'HR');
```

---

### 4️⃣ UPDATE

修改現有資料。

```sql
UPDATE employees SET department = 'Finance' WHERE name = 'Alice';
```

---

### 5️⃣ DELETE

刪除資料列。

```sql
DELETE FROM employees WHERE age < 25;
```

> ⚠️ 若省略 `WHERE`，會刪除整張表的所有資料。

---

### 6️⃣ CREATE TABLE

建立新資料表。

```sql
CREATE TABLE employees (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL,
  age INTEGER,
  department TEXT
);
```

---

### 7️⃣ ALTER TABLE

修改資料表結構。

新增欄位：

```sql
ALTER TABLE employees ADD COLUMN salary INTEGER;
```

---

### 8️⃣ DROP TABLE

刪除整個資料表（包含資料）。

```sql
DROP TABLE employees;
```

---

## 🔍 查詢與篩選語法

### WHERE 條件

可使用邏輯運算子：

```sql
SELECT * FROM employees WHERE department = 'IT' AND age > 25;
```

### ORDER BY 排序

```sql
SELECT * FROM employees ORDER BY age DESC;
```

### LIMIT 限制筆數

```sql
SELECT * FROM employees LIMIT 5;
```

---

## 🔎 模糊查詢 — LIKE

`LIKE` 搭配萬用字元 `%` 或 `_` 進行模糊比對。

| 符號  | 功能       | 範例                                   |
| --- | -------- | ------------------------------------ |
| `%` | 代表任意多個字元 | `WHERE name LIKE 'A%'`（以 A 開頭）       |
| `_` | 代表任意單一字元 | `WHERE name LIKE '_im'`（第二與第三字元為 im） |

例：

```sql
SELECT * FROM employees WHERE name LIKE '%son';  -- 以 son 結尾
```

---

## 🔢 範圍查詢 — BETWEEN

用於查詢介於兩個值之間的資料。

```sql
SELECT * FROM employees WHERE age BETWEEN 25 AND 40;
```

可與日期或文字使用：

```sql
SELECT * FROM orders WHERE order_date BETWEEN '2024-01-01' AND '2024-06-30';
```

---

## 🧩 條件式 — CASE WHEN

`CASE WHEN` 用於在查詢中進行條件判斷。

```sql
SELECT name,
       age,
       CASE
           WHEN age < 30 THEN 'Young'
           WHEN age BETWEEN 30 AND 50 THEN 'Middle-aged'
           ELSE 'Senior'
       END AS age_group
FROM employees;
```

也可以用在聚合運算：

```sql
SELECT department,
       SUM(CASE WHEN gender = 'F' THEN 1 ELSE 0 END) AS female_count
FROM employees
GROUP BY department;
```

---

## 🧮 聚合函數 (Aggregate Functions)

| 函數      | 功能   |
| ------- | ---- |
| COUNT() | 計算筆數 |
| SUM()   | 加總   |
| AVG()   | 平均值  |
| MIN()   | 最小值  |
| MAX()   | 最大值  |

例：

```sql
SELECT department, AVG(salary) FROM employees GROUP BY department;
```

---

## 🔗 表格關聯 (JOIN)

```sql
SELECT employees.name, departments.department_name
FROM employees
JOIN departments ON employees.department_id = departments.id;
```

JOIN 類型：

* **INNER JOIN**：僅回傳兩邊都有對應的資料。
* **LEFT JOIN**：回傳左表所有資料，即使右表無對應值。
* **RIGHT JOIN**：回傳右表所有資料，即使左表無對應值。
* **FULL JOIN**：合併兩表所有資料。

---

## 🧠 總結

SQL 是資料分析與後端開發中最常用的語言之一。熟悉基本語法（`SELECT`、`INSERT`、`UPDATE`、`DELETE`）與進階查詢（`JOIN`、`CASE WHEN`、`LIKE`、`BETWEEN`、`DISTINCT`）即可應付多數實務情境。
