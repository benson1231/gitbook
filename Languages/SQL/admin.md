# 🗄️ PostgreSQL 權限與角色管理筆記

本文整理 PostgreSQL 權限控制（Privileges）與角色管理（Roles）的核心概念與常見操作，包括角色新增/刪除、權限授予/回收，以及行級（Row-level）與列級（Column-level）安全控制。

---

## 🧩 基本概念

### 🔹 使用者與角色（Roles）

PostgreSQL 將「使用者」與「角色」統一成 *Role* 概念。每個角色可：

* 登入（`LOGIN`）
* 建立資料庫（`CREATEDB`）
* 建立其他角色（`CREATEROLE`）
* 擁有超級權限（`SUPERUSER`）

```sql
-- 建立可登入的角色
CREATE ROLE dev_user LOGIN PASSWORD 'mypassword';

-- 建立管理者角色
CREATE ROLE admin LOGIN CREATEDB CREATEROLE;
```

### 🔹 檢視現有角色

```sql
\du  -- psql 指令，顯示所有角色與屬性
```

### 🔹 刪除角色

```sql
DROP ROLE dev_user;
```

> ⚠️ 若角色仍擁有資料表或其他物件，需先轉移所有權再刪除。

```sql
REASSIGN OWNED BY dev_user TO postgres;
DROP OWNED BY dev_user;
DROP ROLE dev_user;
```

---

## 🔐 權限（Privileges）

PostgreSQL 權限可分為三層：

1. **資料庫層級（Database-level）**
2. **資料表與檢視層級（Table/View-level）**
3. **欄位（Column-level）與行級安全（RLS）**

### 🔸 常見權限種類

| 權限        | 適用物件             | 功能                   |
| --------- | ---------------- | -------------------- |
| `CONNECT` | Database         | 允許連線到資料庫             |
| `CREATE`  | Schema, Database | 可建立物件（table, view 等） |
| `SELECT`  | Table, View      | 查詢資料                 |
| `INSERT`  | Table            | 新增資料                 |
| `UPDATE`  | Table, Column    | 修改資料（可限欄位）           |
| `DELETE`  | Table            | 刪除資料                 |
| `USAGE`   | Schema, Sequence | 使用 schema 或序列        |
| `EXECUTE` | Function         | 執行函式                 |

### 🔸 授予權限（GRANT）

```sql
-- 授予使用者查詢與新增權限
GRANT SELECT, INSERT ON employees TO dev_user;

-- 授予整個 schema 權限
GRANT USAGE, CREATE ON SCHEMA public TO dev_user;

-- 授予整個資料庫連線權限
GRANT CONNECT ON DATABASE company_db TO dev_user;
```

### 🔸 回收權限（REVOKE）

```sql
REVOKE INSERT ON employees FROM dev_user;
REVOKE ALL PRIVILEGES ON DATABASE company_db FROM dev_user;
```

---

## 🧱 欄位（Column-level）權限

你可以只允許使用者讀取或修改特定欄位。

```sql
-- 只允許讀取 id 與 name 欄位
GRANT SELECT (id, name) ON employees TO dev_user;

-- 允許修改 email 欄位
GRANT UPDATE (email) ON employees TO dev_user;
```

---

## 🧬 行級安全（Row-Level Security, RLS）

PostgreSQL 支援行級安全策略，讓不同角色只能看到或操作特定的資料列。

### 啟用行級安全

```sql
ALTER TABLE employees ENABLE ROW LEVEL SECURITY;
```

### 建立行級策略（POLICY）

```sql
-- 僅允許使用者看到自己部門的資料
CREATE POLICY same_department_policy
ON employees
FOR SELECT
USING (department_id = current_setting('app.current_department')::int);

-- 限制新增資料的範圍
CREATE POLICY insert_own_department
ON employees
FOR INSERT
WITH CHECK (department_id = current_setting('app.current_department')::int);
```

### 管理策略

```sql
-- 查看策略
\d+ employees  -- 會顯示啟用的 policy 列表

-- 刪除策略
DROP POLICY same_department_policy ON employees;

-- 停用 RLS
ALTER TABLE employees DISABLE ROW LEVEL SECURITY;
```

---

## 🧰 權限綜合範例

假設我們有一個公司資料庫 `company_db`，要建立一位只能看自己部門員工的使用者：

```sql
-- 建立角色
CREATE ROLE hr_user LOGIN PASSWORD '1234';

-- 授予連線與讀取權限
GRANT CONNECT ON DATABASE company_db TO hr_user;
GRANT USAGE ON SCHEMA public TO hr_user;
GRANT SELECT ON employees TO hr_user;

-- 啟用 RLS 並設定策略
ALTER TABLE employees ENABLE ROW LEVEL SECURITY;
CREATE POLICY hr_self_view
ON employees
FOR SELECT
USING (department_id = current_setting('app.current_department')::int);
```

> ✅ 現在 `hr_user` 只會看到與自己部門相符的員工資料。

---

## 🧾 管理技巧

* 使用 `ALTER ROLE` 修改角色屬性：

  ```sql
  ALTER ROLE dev_user WITH PASSWORD 'newpass' LOGIN;
  ALTER ROLE hr_user SET search_path = 'public';
  ```
* 使用 `pg_roles` 或 `pg_user` 系統視圖查詢角色屬性：

  ```sql
  SELECT rolname, rolsuper, rolcreatedb, rolcreaterole, rolcanlogin FROM pg_roles;
  ```
* 使用 `\z`（或 `\dp`）檢查資料表的權限分佈。

---

## 📚 延伸閱讀

* [PostgreSQL 官方文件 – Privileges](https://www.postgresql.org/docs/current/ddl-priv.html)
* [PostgreSQL 官方文件 – Row Security Policies](https://www.postgresql.org/docs/current/ddl-rowsecurity.html)
* [psql 命令快速查詢權限與角色](https://www.postgresql.org/docs/current/app-psql.html)

---

💡 **建議**：在開發或生產環境，將權限分層為 `readonly`、`readwrite`、`admin` 等角色模板，避免直接給最終使用者廣泛權限，以降低資料洩漏與誤操作風險。
