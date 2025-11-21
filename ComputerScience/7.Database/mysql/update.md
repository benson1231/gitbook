# MySQL 更新與刪除資料教學

本篇介紹如何使用 `UPDATE` 與 `DELETE` 指令來修改與刪除資料。

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

## 更新資料（`UPDATE`）

### 更新單筆資料

```sql
-- 將 id 為 1 的使用者年齡改為 21
UPDATE mytable SET age = 21 WHERE id = 1;
```

### 同時更新多個欄位

```sql
-- 更新 name 與 age
UPDATE mytable SET name = 'benson_wu', age = 22 WHERE id = 1;
```

### 沒有加 WHERE（會更新所有資料）

```sql
-- 小心！以下語法會修改整個資料表的 age 欄位
UPDATE mytable SET age = 99;
```

## 刪除資料（`DELETE`）

### 刪除指定資料

```sql
-- 刪除 name 為 'alice' 的紀錄
DELETE FROM mytable WHERE name = 'alice';
```

### 沒有加 WHERE（會刪除所有資料）

```sql
-- 小心！以下語法會清空整張資料表
DELETE FROM mytable;
```

---

建議在執行 `UPDATE` 或 `DELETE` 指令前先用 `SELECT` 確認條件正確，避免誤修改或刪除資料。
