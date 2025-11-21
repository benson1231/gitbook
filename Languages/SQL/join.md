# SQL JOIN 與 UNION 教學

## 🧩 JOIN：資料表關聯

`JOIN` 用於將多個資料表中的資料以關聯鍵（key）結合起來。

### 🔹 INNER JOIN

回傳兩個資料表中符合條件的資料。

```sql
SELECT employees.name, departments.department_name
FROM employees
INNER JOIN departments
ON employees.department_id = departments.id;
```

➡️ 只顯示同時存在於 `employees` 與 `departments` 的匹配資料。

---

### 🔹 LEFT JOIN

回傳左表的所有資料，即使右表沒有匹配值。

```sql
SELECT employees.name, departments.department_name
FROM employees
LEFT JOIN departments
ON employees.department_id = departments.id;
```

➡️ 若右表無對應資料，該欄位以 `NULL` 顯示。

---

### 🔹 RIGHT JOIN

回傳右表的所有資料，即使左表沒有匹配值。

```sql
SELECT employees.name, departments.department_name
FROM employees
RIGHT JOIN departments
ON employees.department_id = departments.id;
```

---

### 🔹 FULL JOIN

合併兩表所有資料，若一方無對應則以 `NULL` 補上。

```sql
SELECT employees.name, departments.department_name
FROM employees
FULL JOIN departments
ON employees.department_id = departments.id;
```

---

### 🔹 CROSS JOIN

產生兩表的笛卡兒積（每筆左表資料對應右表所有資料）。

```sql
SELECT * FROM employees CROSS JOIN departments;
```

➡️ 若左表有 10 筆、右表有 5 筆，結果為 50 筆。

---

### 🔸 JOIN 類型比較

| JOIN 類型    | 回傳內容        | NULL 顯示     |
| ---------- | ----------- | ----------- |
| INNER JOIN | 兩表都有的資料     | 無           |
| LEFT JOIN  | 左表全部 + 右表匹配 | 右表缺值為 NULL  |
| RIGHT JOIN | 右表全部 + 左表匹配 | 左表缺值為 NULL  |
| FULL JOIN  | 合併左右兩表全部    | 缺值以 NULL 顯示 |
| CROSS JOIN | 左表 × 右表組合   | 無條件匹配       |

---

## 🔗 UNION：合併查詢結果

### 🔹 UNION

將兩個查詢結果合併，並**移除重複值**。

```sql
SELECT name FROM employees
UNION
SELECT name FROM managers;
```

➡️ 結果中相同的名稱僅會出現一次。

---

### 🔹 UNION ALL

將結果合併但**保留重複值**。

```sql
SELECT name FROM employees
UNION ALL
SELECT name FROM managers;
```

---

### 🔸 UNION 與 UNION ALL 比較

| 指令        | 是否去除重複 | 效能      |
| --------- | ------ | ------- |
| UNION     | ✅ 是    | 較慢（需去重） |
| UNION ALL | ❌ 否    | 較快      |

---

## 🧱 WITH（CTE，共用表表達式）

`WITH` 用於建立臨時查詢結果，可被後續 SQL 使用，使查詢更清晰。

### 🔹 基本語法

```sql
WITH department_salary AS (
  SELECT department_id, AVG(salary) AS avg_salary
  FROM employees
  GROUP BY department_id
)
SELECT e.name, d.avg_salary
FROM employees e
JOIN department_salary d
ON e.department_id = d.department_id;
```

➡️ 此語法先建立名為 `department_salary` 的臨時結果，再於主查詢中重複使用。

---

## 🧠 小結

* **JOIN** 用於橫向合併資料表（根據欄位關聯）。
* **UNION / UNION ALL** 用於縱向合併查詢結果。
* **WITH (CTE)** 可重用中間查詢，讓 SQL 更具可讀性。

👉 建議在實務應用中搭配使用 JOIN、UNION 與 CTE，以組合出結構化且高效的資料查詢流程。
