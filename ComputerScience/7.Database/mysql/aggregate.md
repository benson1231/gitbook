# MySQL 聚合函式（Aggregate Functions）教學

聚合函式用於統計資料表中的多筆資料，常搭配 `GROUP BY` 使用，適合進行彙總與統計查詢。

## 常見聚合函式

| 函式名稱                 | 說明                        |
| -------------------- | ------------------------- |
| `COUNT()`            | 計算筆數                      |
| `SUM()`              | 加總欄位值                     |
| `AVG()`              | 計算平均值                     |
| `MIN()`              | 取得最小值                     |
| `MAX()`              | 取得最大值                     |
| `STD()` / `STDDEV()` | 計算標準差（Standard Deviation） |

## 基本語法範例

假設有一張 `orders` 資料表：

```sql
CREATE TABLE orders (
    id INT PRIMARY KEY,
    customer VARCHAR(50),
    amount DECIMAL(10,2)
);

INSERT INTO orders VALUES
(1, 'Alice', 120.50),
(2, 'Bob', 89.00),
(3, 'Alice', 55.00),
(4, 'Charlie', 140.75);
```

### 1. 總筆數

```sql
SELECT COUNT(*) FROM orders;
```

### 2. 金額總和

```sql
SELECT SUM(amount) FROM orders;
```

### 3. 平均金額

```sql
SELECT AVG(amount) FROM orders;
```

### 4. 最大與最小金額

```sql
SELECT MAX(amount), MIN(amount) FROM orders;
```

### 5. 標準差

```sql
SELECT STD(amount) FROM orders;
-- 或使用 STDDEV(amount)
```

## 搭配 GROUP BY

```sql
-- 各客戶的總金額、訂單數量與金額標準差
SELECT customer,
       COUNT(*) AS order_count,
       SUM(amount) AS total_amount,
       STD(amount) AS std_amount
FROM orders
GROUP BY customer;
```

## 搭配 HAVING 篩選彙總結果

```sql
-- 查詢消費總額超過 100 的客戶
SELECT customer, SUM(amount) AS total_amount
FROM orders
GROUP BY customer
HAVING total_amount > 100;
```

---

聚合函式是商業報表與資料分析中不可或缺的工具，善用 `GROUP BY` 與 `HAVING` 能大幅強化查詢彈性與表達力。

> 📌 備註：`STD()` 與 `STDDEV()` 在 MySQL 中功能相同，用於計算樣本標準差，適用於統計分析。
