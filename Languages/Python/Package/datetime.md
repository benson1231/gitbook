# 🕒 Python `datetime` 模組教學

Python 的 `datetime` 模組提供了日期與時間的操作功能，能用於時間計算、格式化、時區轉換等，是資料分析、排程與時間標註常用的工具。

---

## 📦 匯入模組

```python
import datetime
```

---

## 🧩 主要類別

`datetime` 模組包含幾個核心類別：

| 類別          | 說明              |
| ----------- | --------------- |
| `date`      | 表示日期（年、月、日）     |
| `time`      | 表示時間（時、分、秒、微秒）  |
| `datetime`  | 同時包含日期與時間       |
| `timedelta` | 表示時間差，用於日期/時間運算 |
| `timezone`  | 表示時區資訊          |

---

## 📅 建立日期與時間物件

### 1️⃣ 建立日期 (`date`)

```python
from datetime import date

# 建立特定日期
my_date = date(2025, 10, 8)
print(my_date)  # 2025-10-08

# 取得今天日期
today = date.today()
print(today)
```

### 2️⃣ 建立時間 (`time`)

```python
from datetime import time

my_time = time(14, 30, 15)
print(my_time)  # 14:30:15
```

### 3️⃣ 建立日期時間 (`datetime`)

```python
from datetime import datetime

now = datetime.now()
print(now)  # 2025-10-08 14:30:15.123456

specific_dt = datetime(2025, 10, 8, 14, 30, 15)
print(specific_dt)
```

---

## 🔢 日期與時間的屬性

```python
dt = datetime(2025, 10, 8, 14, 30, 15)
print(dt.year)   # 2025
print(dt.month)  # 10
print(dt.day)    # 8
print(dt.hour)   # 14
print(dt.minute) # 30
```

---

## 🧮 日期與時間運算 (`timedelta`)

```python
from datetime import timedelta

now = datetime.now()
one_week = timedelta(days=7)

print(now + one_week)   # 一週後的日期時間
print(now - one_week)   # 一週前的日期時間
```

也可用於時間差計算：

```python
d1 = datetime(2025, 10, 8)
d2 = datetime(2025, 9, 1)

diff = d1 - d2
print(diff.days)  # 37
```

---

## 🧭 時間格式化 (`strftime` / `strptime`)

### 格式化輸出 (`strftime`)

```python
now = datetime.now()
print(now.strftime("%Y-%m-%d %H:%M:%S"))  # 2025-10-08 14:30:15
```

### 解析字串成時間 (`strptime`)

```python
dt_str = "2025-10-08 14:30:15"
dt = datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
print(dt)
```

### 常用格式符號

| 格式   | 意義        | 範例        |
| ---- | --------- | --------- |
| `%Y` | 西元年 (四位數) | 2025      |
| `%m` | 月 (01–12) | 10        |
| `%d` | 日 (01–31) | 08        |
| `%H` | 時 (00–23) | 14        |
| `%M` | 分 (00–59) | 30        |
| `%S` | 秒 (00–59) | 15        |
| `%A` | 星期名稱      | Wednesday |
| `%a` | 星期縮寫      | Wed       |
| `%B` | 月份名稱      | October   |

---

## 🌍 時區 (`timezone`)

```python
from datetime import datetime, timezone, timedelta

utc_now = datetime.now(timezone.utc)
print(utc_now)  # UTC 時間

taipei_tz = timezone(timedelta(hours=8))
taipei_now = datetime.now(taipei_tz)
print(taipei_now)
```

---

## 💡 小技巧

1. 比較時間時，務必在相同時區下比較。
2. 若僅需日期部分，可用 `.date()` 取得。
3. 可使用 `isoformat()` 與 `fromisoformat()` 進行 ISO 標準轉換。

   ```python
   iso_str = now.isoformat()
   print(iso_str)
   parsed = datetime.fromisoformat(iso_str)
   ```

---

## ✅ 總結

`datetime` 模組整合了日期、時間、時區與格式化操作，是處理時間資料的核心工具。掌握 `datetime`、`timedelta`、`strftime`、`strptime` 這幾個關鍵功能，就能輕鬆應付多數時間處理任務。
