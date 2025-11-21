# 使用 mysql-connector-python 操作 MySQL 資料庫

`mysql-connector-python` 是由 Oracle 提供的官方 Python 套件，用於與 MySQL 資料庫進行連線與操作。

## 安裝

```bash
pip install mysql-connector-python
```

## 基本連線與查詢範例

```python
import mysql.connector

# 建立連線
conn = mysql.connector.connect(
    host="localhost",
    user="root",
    password="your_password",
    database="test"
)

cursor = conn.cursor()

# 執行查詢
cursor.execute("SELECT * FROM mytable")

# 取得結果
results = cursor.fetchall()
for row in results:
    print(row)

# 關閉連線
cursor.close()
conn.close()
```

## 建立資料表與新增資料

```python
cursor.execute("""
CREATE TABLE IF NOT EXISTS mytable (
    id INT AUTO_INCREMENT PRIMARY KEY,
    name VARCHAR(100),
    age INT
)
""")

cursor.execute("INSERT INTO mytable (name, age) VALUES (%s, %s)", ("Alice", 30))
conn.commit()
```

## 使用參數化查詢（防止 SQL Injection）

```python
name = "Bob"
age = 25
cursor.execute("INSERT INTO mytable (name, age) VALUES (%s, %s)", (name, age))
conn.commit()
```

## 查詢與條件篩選

```python
cursor.execute("SELECT * FROM mytable WHERE age > %s", (20,))
for row in cursor.fetchall():
    print(row)
```

## 更新與刪除資料

```python
# 更新資料
cursor.execute("UPDATE mytable SET age = %s WHERE name = %s", (28, "Alice"))
conn.commit()

# 刪除資料
cursor.execute("DELETE FROM mytable WHERE name = %s", ("Bob",))
conn.commit()
```

## 取得資料筆數與逐筆讀取

```python
cursor.execute("SELECT * FROM mytable")

# 總筆數
print("rows:", cursor.rowcount)

# 逐筆讀取
row = cursor.fetchone()
while row:
    print(row)
    row = cursor.fetchone()
```

```python
cursor.execute("SELECT * FROM mytable WHERE age > %s", (20,))
for row in cursor.fetchall():
print(row)
```

## 錯誤處理範例
```python
try:
    conn = mysql.connector.connect(...)
    cursor = conn.cursor()
    # 執行動作
except mysql.connector.Error as err:
    print(f"MySQL 錯誤: {err}")
finally:
    if cursor:
        cursor.close()
    if conn:
        conn.close()
```

---

`mysql-connector-python` 適合用於資料存取、後端整合與自動化資料處理。也可搭配 `pandas` 使用 `pd.read_sql()` 整合分析流程。
