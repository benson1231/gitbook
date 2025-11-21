# SQL 聚合函數（Aggregate Functions）教學

## 🧭 簡介

**聚合函數 (Aggregate Functions)** 是 SQL 中用來對多筆資料進行運算、彙總或統計的函數。常見用途包括計算筆數、平均值、最大值、最小值與總和。

聚合函數通常與 `GROUP BY` 子句搭配使用，用於依群組彙總資料，例如計算各部門的平均薪資或每月銷售總額。

---

## 📘 常見聚合函數

| 函數名稱      | 功能說明 | 範例                                   |
| --------- | ---- | ------------------------------------ |
| `COUNT()` | 計算筆數 | `SELECT COUNT(*) FROM employees;`    |
| `SUM()`   | 加總   | `SELECT SUM(salary) FROM employees;` |
| `AVG()`   | 平均值  | `SELECT AVG(salary) FROM employees;` |
| `MIN()`   | 最小值  | `SELECT MIN(age) FROM employees;`    |
| `MAX()`   | 最大值  | `SELECT MAX(age) FROM employees;`    |

---

## 📊 COUNT()

`COUNT()` 用於計算資料筆數，可搭配 `*` 或指定欄位使用。

```sql
-- 計算所有員工數
SELECT COUNT(*) FROM employees;

-- 計算指定部門的員工數
SELECT COUNT(*) FROM employees WHERE department = 'IT';
```

> 💡 `COUNT(column_name)` 不會計算 NULL 值。

---

## 💰 SUM()

`SUM()` 用於數值欄位的總和。

```sql
-- 計算所有員工的總薪資
SELECT SUM(salary) FROM employees;

-- 計算每個部門的總薪資
SELECT department, SUM(salary) FROM employees GROUP BY department;
```

---

## 📈 AVG()

`AVG()` 計算平均值。

```sql
-- 全體員工平均薪資
SELECT AVG(salary) FROM employees;

-- 各部門平均薪資
SELECT department, AVG(salary) FROM employees GROUP BY department;
```

---

## 📉 MIN() 與 MAX()

`MIN()` 取得最小值，`MAX()` 取得最大值。

```sql
-- 最年輕與最年長員工年齡
SELECT MIN(age) AS youngest, MAX(age) AS oldest FROM employees;
```

---

## 🧩 GROUP BY 群組聚合

`GROUP BY` 用於依照一個或多個欄位分組，再進行彙總。

```sql
-- 各部門的員工數與平均薪資
SELECT department, COUNT(*) AS total_employees, AVG(salary) AS avg_salary
FROM employees
GROUP BY department;
```

> ⚠️ `GROUP BY` 後的 SELECT 子句中，只能包含群組欄位與聚合函數。

---

## ⚙️ HAVING 條件過濾群組

`HAVING` 用於篩選群組結果（不同於 `WHERE` 是篩選個別資料列）。

```sql
-- 顯示平均薪資超過 50000 的部門
SELECT department, AVG(salary) AS avg_salary
FROM employees
GROUP BY department
HAVING AVG(salary) > 50000;
```

---

## 🧠 多層聚合與運算

聚合結果可以再進行數學運算或巢狀查詢：

```sql
-- 計算所有部門平均薪資的平均值
SELECT AVG(avg_salary)
FROM (
  SELECT AVG(salary) AS avg_salary
  FROM employees
  GROUP BY department
) AS dept_avg;
```

---

## 🎯 小結

* 聚合函數能快速進行資料統計與彙總。
* 常與 `GROUP BY`、`HAVING` 一起使用。
* 聚合結果可搭配別名（`AS`）提升可讀性。
* 聚合函數會忽略 `NULL` 值（除非使用 `COUNT(*)`）。

---

### ✅ 範例練習

```sql
-- 計算各部門的員工人數與最高薪資
SELECT department, COUNT(*) AS num_employees, MAX(salary) AS top_salary
FROM employees
GROUP BY department
ORDER BY top_salary DESC;
```
