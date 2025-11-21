# MySQL 外鍵（FOREIGN KEY）教學

外鍵用來建立資料表之間的參照完整性（Referential Integrity），確保關聯欄位的值在父表（Parent Table）中存在。

## 基本概念

* **父表 (Parent)**：提供主鍵，被其他表參照。
* **子表 (Child)**：包含外鍵欄位，參照父表的主鍵。
* **參照行為 (Referential Actions)**：定義父表資料被更新或刪除時，子表資料的處理方式，例如 `CASCADE`、`SET NULL`、`RESTRICT`。

## 資料表範例（建立時即指定外鍵）

```sql
-- 父表：departments
CREATE TABLE departments (
    dept_id INT PRIMARY KEY,
    dept_name VARCHAR(100)
);

-- 子表：employees
CREATE TABLE employees (
    emp_id INT PRIMARY KEY,
    emp_name VARCHAR(100),
    dept_id INT,
    CONSTRAINT fk_dept
        FOREIGN KEY (dept_id)
        REFERENCES departments(dept_id)
        ON DELETE CASCADE
        ON UPDATE CASCADE
);
```

### 選項說明

* `ON DELETE CASCADE`：當父表紀錄被刪除時，自動刪除子表中對應紀錄。
* `ON UPDATE CASCADE`：當父表主鍵更新時，子表外鍵同步更新。
* 其他選項：`SET NULL`、`RESTRICT`、`NO ACTION`。

## 新增外鍵到既有資料表

```sql
ALTER TABLE employees
ADD CONSTRAINT fk_dept
FOREIGN KEY (dept_id)
REFERENCES departments(dept_id)
ON DELETE SET NULL;
```

## 移除外鍵

```sql
ALTER TABLE employees
DROP FOREIGN KEY fk_dept;
```

## 注意事項

1. 子表外鍵欄位必須有索引（InnoDB 會自動建立）。
2. 外鍵與參照欄位的資料型別與長度須相同。
3. 僅支援 InnoDB 引擎；MyISAM 不支援外鍵。

---

外鍵可保護資料一致性，但亦可能影響批次刪除與更新效能。如需大量操作，請先檢視 `FOREIGN_KEY_CHECKS` 或採用事務控制 (`START TRANSACTION`)。
