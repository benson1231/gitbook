# 🧱 SQL 資料型別（Data Types）教學

SQL 資料型別（Data Types）用來定義每個欄位可儲存的資料種類，例如整數、文字、日期或浮點數。不同資料庫（如 MySQL、PostgreSQL、SQLite）在型別名稱上略有不同，但概念大致相同。

---

## 🧩 數值型（Numeric Types）

| 型別                                | 說明                      | 範例                            |
| --------------------------------- | ----------------------- | ----------------------------- |
| `INT` / `INTEGER`                 | 整數（無小數）                 | `100`, `-25`                  |
| `SMALLINT`                        | 小範圍整數                   | `32767`, `-32768`             |
| `BIGINT`                          | 大範圍整數                   | `9223372036854775807`         |
| `DECIMAL(p, s)` / `NUMERIC(p, s)` | 精確小數，`p` 為總位數，`s` 為小數位數 | `DECIMAL(10, 2)` → `12345.67` |
| `FLOAT` / `REAL` / `DOUBLE`       | 浮點數（近似值）                | `3.14159`, `-0.001`           |

📘 **建議使用**：金錢或需要高精度的數值使用 `DECIMAL`，計算用數值可用 `FLOAT`。

---

## 🔤 文字型（Character Types）

| 型別           | 說明               | 範例                              |
| ------------ | ---------------- | ------------------------------- |
| `CHAR(n)`    | 固定長度字串，不足部分以空白補齊 | `'YES '` (CHAR(4))              |
| `VARCHAR(n)` | 可變長度字串           | `'HELLO'` (VARCHAR(10))         |
| `TEXT`       | 長文字資料            | `'This is a long paragraph...'` |

📘 **建議使用**：

* 固定長度欄位（如性別代碼）用 `CHAR`。
* 一般字串（如姓名、Email）用 `VARCHAR`。
* 內容較長（如備註、描述）用 `TEXT`。

---

## 📅 日期與時間型（Date & Time Types）

| 型別                       | 說明             | 範例                      |
| ------------------------ | -------------- | ----------------------- |
| `DATE`                   | 日期（YYYY-MM-DD） | `'2025-10-08'`          |
| `TIME`                   | 時間（HH:MI:SS）   | `'14:30:00'`            |
| `DATETIME` / `TIMESTAMP` | 日期 + 時間        | `'2025-10-08 14:30:00'` |
| `YEAR`                   | 年份             | `'2025'`                |

📘 **注意**：`TIMESTAMP` 通常包含自動更新功能（記錄修改時間），適用於審計紀錄。

---

## 🧮 布林值型（Boolean Type）

| 型別        | 說明                | 範例              |
| --------- | ----------------- | --------------- |
| `BOOLEAN` | 布林值（True / False） | `TRUE`, `FALSE` |

在某些資料庫中（如 MySQL），`BOOLEAN` 實際上是 `TINYINT(1)`，即 `1` 表示真，`0` 表示假。

---

## 📦 其他特殊型別

| 型別     | 說明                                 | 範例                                       |
| ------ | ---------------------------------- | ---------------------------------------- |
| `BLOB` | 二進位資料（Binary Large Object），儲存圖片或檔案 | 影像、音訊檔                                   |
| `JSON` | 結構化 JSON 資料（PostgreSQL / MySQL 支援） | `'[{"id": 1, "name": "Alice"}]'`         |
| `UUID` | 通用唯一識別碼                            | `'550e8400-e29b-41d4-a716-446655440000'` |

---

## 🧠 範例：建立資料表時指定型別

```sql
CREATE TABLE employees (
  id SERIAL PRIMARY KEY,
  name VARCHAR(100) NOT NULL,
  age INT CHECK (age > 0),
  salary DECIMAL(10, 2) DEFAULT 0.00,
  is_active BOOLEAN DEFAULT TRUE,
  hire_date DATE NOT NULL,
  profile JSON
);
```

---

## 🧾 小結

* 數值型（`INT`, `DECIMAL`, `FLOAT`）用於數字資料。
* 文字型（`CHAR`, `VARCHAR`, `TEXT`）用於字串。
* 日期型（`DATE`, `TIME`, `TIMESTAMP`）用於時間資料。
* 布林型（`BOOLEAN`）用於邏輯值。
* 特殊型（`BLOB`, `JSON`, `UUID`）用於結構化或二進位資料。

📘 **選型原則**：

* 根據資料特性與精度需求選擇適當型別。
* 避免使用過大型別以節省儲存空間。
* 在 PostgreSQL、MySQL、SQLite 間遷移時，確認型別支援差異。
