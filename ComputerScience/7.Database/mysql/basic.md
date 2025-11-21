# MySQL basic tutorial

以下為 MySQL 資料庫的基本操作指令，適用於初學者練習建立資料庫、資料表以及進行資料新增與查詢。

### 顯示所有資料庫

```sql
SHOW DATABASES;
```

### 建立資料庫 `test`

```sql
CREATE DATABASE test;
```

### 刪除資料庫 `test`

```sql
DROP DATABASE test;
```

### 使用資料庫 `test`

```sql
USE test;
```

### 顯示目前資料庫中的所有資料表

```sql
SHOW TABLES;
```

### 建立資料表 `mytable`（基本版）

```sql
CREATE TABLE mytable(
    age INT,
    name VARCHAR(255)
);
```

### 建立資料表 `mytable`（進階版）

包含主鍵、自動遞增欄位、NOT NULL 與預設值設定：

```sql
CREATE TABLE mytable(
    id INT PRIMARY KEY AUTO_INCREMENT,
    age INT NOT NULL DEFAULT 18,
    name VARCHAR(255) NOT NULL
);
```

### 插入資料
```sql
INSERT INTO mytable(age, name) VALUES (18, "benson");
INSERT INTO mytable(age, name) VALUES (18, "evan");
```

### 查詢資料

```sql
SELECT * FROM mytable;
```

## 刪除資料表

```sql
DROP TABLE mytable;
```
