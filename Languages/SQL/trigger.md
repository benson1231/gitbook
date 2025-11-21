# 🧩 SQL 觸發器 (Trigger) 教學

## 🧭 什麼是觸發器？

**觸發器（Trigger）** 是一種特殊的資料庫物件，當某個事件（例如 `INSERT`、`UPDATE` 或 `DELETE`）在指定的資料表上發生時，會自動執行預先定義的 SQL 程式。

用途包括：

* 自動記錄資料變更（例如紀錄修改者、修改時間）
* 自動同步資料（例如更新其他表）
* 資料驗證（阻止不合法更新）

---

## ⚙️ 基本語法

```sql
CREATE TRIGGER trigger_name
{BEFORE | AFTER} {INSERT | UPDATE | DELETE}
ON table_name
FOR EACH ROW
BEGIN
    -- 觸發時執行的 SQL 敘述
END;
```

### 🔸 說明：

| 關鍵字            | 意義                            |
| -------------- | ----------------------------- |
| `BEFORE`       | 在事件發生之前執行觸發器                  |
| `AFTER`        | 在事件發生之後執行觸發器                  |
| `FOR EACH ROW` | 表示針對被影響的每一列都執行觸發動作            |
| `NEW`          | 代表新資料列（`INSERT` 或 `UPDATE` 後） |
| `OLD`          | 代表舊資料列（`UPDATE` 或 `DELETE` 前） |

---

## 🧱 實際範例

### 1️⃣ 自動記錄更新時間

```sql
CREATE TRIGGER update_timestamp
BEFORE UPDATE ON employees
FOR EACH ROW
BEGIN
    SET NEW.updated_at = CURRENT_TIMESTAMP;
END;
```

> 🔹 功能：在每次更新員工資料時，自動更新 `updated_at` 欄位。

---

### 2️⃣ 防止刪除重要資料

```sql
CREATE TRIGGER prevent_admin_delete
BEFORE DELETE ON users
FOR EACH ROW
BEGIN
    IF OLD.role = 'admin' THEN
        SIGNAL SQLSTATE '45000'
        SET MESSAGE_TEXT = 'Cannot delete admin accounts!';
    END IF;
END;
```

> 🔹 功能：阻止刪除角色為 `admin` 的使用者。

---

### 3️⃣ 自動同步資料至紀錄表

```sql
CREATE TRIGGER log_employee_insert
AFTER INSERT ON employees
FOR EACH ROW
BEGIN
    INSERT INTO employee_logs (emp_id, action, log_time)
    VALUES (NEW.id, 'INSERT', NOW());
END;
```

> 🔹 功能：每當新增員工時，自動將紀錄寫入 `employee_logs` 表。

---

## 🔍 查看與刪除觸發器

### 查看所有觸發器

```sql
SHOW TRIGGERS;
```

### 刪除觸發器

```sql
DROP TRIGGER IF EXISTS trigger_name;
```

---

## 🧠 小技巧與注意事項

* `BEFORE` 觸發器可修改 `NEW` 值，`AFTER` 則不可。
* 一個資料表可有多個觸發器，但執行順序需明確設計。
* 避免在觸發器中執行過多複雜邏輯，否則會降低效能。
* 在 PostgreSQL 中，觸發器通常搭配 **觸發函數 (Trigger Function)** 使用。

---

## 🧩 PostgreSQL 範例：觸發函數 + 觸發器

```sql
-- 建立觸發函數
CREATE OR REPLACE FUNCTION log_changes()
RETURNS TRIGGER AS $$
BEGIN
    INSERT INTO audit_log(table_name, operation, changed_at)
    VALUES (TG_TABLE_NAME, TG_OP, NOW());
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- 建立觸發器
CREATE TRIGGER after_update_audit
AFTER UPDATE ON employees
FOR EACH ROW
EXECUTE FUNCTION log_changes();
```

> 🔹 PostgreSQL 使用 `EXECUTE FUNCTION` 來呼叫觸發函數。

---

## 🧾 總結

| 功能    | 關鍵字                          | 範例        |
| ----- | ---------------------------- | --------- |
| 新增觸發器 | `CREATE TRIGGER`             | 自動更新欄位    |
| 刪除觸發器 | `DROP TRIGGER`               | 刪除舊觸發器    |
| 查看觸發器 | `SHOW TRIGGERS`              | 查詢目前所有觸發器 |
| 關聯事件  | `INSERT`, `UPDATE`, `DELETE` | 定義觸發條件    |

觸發器是自動化資料處理的強大工具，能協助確保資料一致性、審核追蹤與邏輯自動化。
