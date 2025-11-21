# String Format

在 Python 中，我們常用 `format()` 方法或 f-string 來插入變數到字串中。這兩種方式都可以讓你以可讀、可維護的方式建立包含變數的文字。

---

## 一、`str.format()` 方法

```python
print("Creating circle with diameter {d}".format(d=diameter))
```

* `{d}` 是一個佔位符，會被 `format()` 中提供的變數取代。
* 你也可以使用索引或位置引數：

```python
print("{0} + {1} = {2}".format(2, 3, 5))
```

---

## 二、f-string（Python 3.6+ 推薦用法）

```python
diameter = 10
print(f"Creating circle with diameter {diameter}")
```

* 將變數直接包在 `{}` 中，前方加上 `f` 表示為 f-string。
* 可直接進行運算與格式設定：

```python
radius = 5
area = 3.14159 * radius ** 2
print(f"Area: {area:.2f}")  # 顯示到小數點第 2 位
```

---

## 三、格式控制符號（格式化標記）

| 語法範例    | 說明                    |
| ------- | --------------------- |
| `:.2f`  | 小數點後 2 位（fixed-point） |
| `:,.0f` | 加入千分位逗號、0 小數          |
| `:>10`  | 靠右對齊，寬度 10            |
| `:<10`  | 靠左對齊，寬度 10            |
| `:^10`  | 置中對齊，寬度 10            |

### 範例：

```python
value = 1234.56789
print(f"{value:.2f}")      # 輸出：1234.57
print(f"{value:,.1f}")     # 輸出：1,234.6
```

---

這些格式化技巧在建立報表、訊息提示、數據輸出時非常實用，能幫助你更清楚地控制文字與數值的輸出樣式。
