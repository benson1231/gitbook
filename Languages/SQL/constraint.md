# SQL Constraint 教學

## 🧭 簡介

**Constraint（約束條件）** 是 SQL 中用來限制資料表中欄位可接受值的規則，用以確保資料的**完整性（integrity）**與**一致性（consistency）**。

常見的 Constraint 包括：

| Constraint  | 功能說明                    |
| ----------- | ----------------------- |
| PRIMARY KEY | 唯一標識每一筆資料，不可重複、不可為 NULL |
| FOREIGN KEY | 建立兩個表格之間的關聯關係           |
| UNIQUE      | 欄位值不可重複，但可為 NULL        |
| NOT NULL    | 欄位不可為 NULL（必須有值）        |
| CHECK       | 限制欄位值必須符合指定條件           |
| DEFAULT     | 指定欄位的預設值                |

---

## 🔑 PRIMARY KEY

唯一標識資料列，每個表只能有一個主鍵。

```sql
CREATE TABLE employees (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL,
  department TEXT
);
```

主鍵也可以設定在多個欄位上（複合主鍵）：

```sql
CREATE TABLE orders (
  order_id INTEGER,
  product_id INTEGER,
  PRIMARY KEY (order_id, product_id)
);
```

---

## 🔗 FOREIGN KEY

用來建立兩張表格間的**參照關係**（referential integrity）。

```sql
CREATE TABLE departments (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL
);

CREATE TABLE employees (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL,
  department_id INTEGER,
  FOREIGN KEY (department_id) REFERENCES departments(id)
);
```

如果刪除或修改被參照的紀錄，可設定相應動作：

```sql
FOREIGN KEY (department_id) REFERENCES departments(id)
  ON DELETE CASCADE
  ON UPDATE SET NULL;
```

| 動作                   | 說明              |
| -------------------- | --------------- |
| CASCADE              | 同步刪除或更新子表資料     |
| SET NULL             | 被參照資料刪除時設為 NULL |
| NO ACTION / RESTRICT | 拒絕操作，維持完整性      |

---

## 🚫 UNIQUE

確保欄位值不重複。

```sql
CREATE TABLE users (
  id INTEGER PRIMARY KEY,
  email TEXT UNIQUE
);
```

可在多個欄位組合上使用：

```sql
CREATE TABLE enrollments (
  student_id INTEGER,
  course_id INTEGER,
  UNIQUE (student_id, course_id)
);
```

---

## ⚠️ NOT NULL

要求欄位必須有值。

```sql
CREATE TABLE employees (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL
);
```

---

## 🧩 CHECK

限制欄位必須符合條件。

```sql
CREATE TABLE employees (
  id INTEGER PRIMARY KEY,
  age INTEGER CHECK (age >= 18)
);
```

也可同時檢查多個欄位：

```sql
CREATE TABLE products (
  id INTEGER PRIMARY KEY,
  price DECIMAL(10,2),
  discount DECIMAL(10,2),
  CHECK (discount <= price)
);
```

---

## 🎯 DEFAULT

指定欄位的預設值。

```sql
CREATE TABLE employees (
  id INTEGER PRIMARY KEY,
  name TEXT NOT NULL,
  status TEXT DEFAULT 'active'
);
```

插入資料時若未指定該欄位，會自動使用預設值。

```sql
INSERT INTO employees (name) VALUES ('Alice');
-- status 會自動設定為 'active'
```

---

## 🧱 ALTER TABLE 新增或移除 Constraint

### 新增 Constraint

```sql
ALTER TABLE employees ADD CONSTRAINT chk_age CHECK (age >= 18);
```

### 刪除 Constraint

不同資料庫系統語法略有不同：

```sql
ALTER TABLE employees DROP CONSTRAINT chk_age;  -- PostgreSQL
ALTER TABLE employees DROP CHECK chk_age;       -- SQLite
```

---

## 🧠 總結

Constraint 是維護資料正確性的重要機制。合理運用 PRIMARY KEY、FOREIGN KEY、UNIQUE、NOT NULL、CHECK 與 DEFAULT，可大幅減少資料錯誤與冗餘，並強化資料庫的可靠性。
