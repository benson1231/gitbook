# 🐍 Python List 教學筆記

## 📘 基本介紹

`list` 是 Python 中最常用的資料結構之一，用於**儲存有序、可變（mutable）**的元素集合。列表可以包含任意型別的資料，例如整數、字串、甚至其他列表。

```python
# 建立一個 list
fruits = ["apple", "banana", "cherry"]

# list 可以混合不同型別
data = [1, "hello", 3.14, True]
```

### 🔹 特性

* **有序（Ordered）**：元素的順序固定，可透過索引（index）訪問。
* **可變（Mutable）**：可以修改、增加或刪除元素。
* **可包含重複元素**。
* **支援巢狀結構（Nested List）**。

---

## 📍 建立與索引操作

```python
# 建立空列表
empty_list = []

# 以 list() 函式建立
nums = list((1, 2, 3))

# 透過索引存取
print(fruits[0])   # apple
print(fruits[-1])  # cherry（反向索引）

# 修改元素
fruits[1] = "mango"

# 切片 (slicing)
print(fruits[0:2])  # ['apple', 'mango']
```

---

## ⚙️ 常用方法與操作

### 1️⃣ 新增元素

```python
# 在尾端新增
fruits.append("orange")

# 插入到指定位置
fruits.insert(1, "grape")

# 合併兩個 list
more = ["kiwi", "melon"]
fruits.extend(more)
```

### 2️⃣ 刪除元素

```python
# 移除指定值
fruits.remove("banana")  # 若不存在會引發 ValueError

# 移除特定位置
fruits.pop(2)  # 回傳並刪除該元素

# 清空整個 list
fruits.clear()
```

### 3️⃣ 搜尋與計算

```python
nums = [1, 2, 3, 2, 4]

print(nums.index(2))  # 找出第一次出現的位置
print(nums.count(2))  # 計算出現次數
```

### 4️⃣ 排序與反轉

```python
nums.sort()          # 就地排序（由小到大）
nums.sort(reverse=True)  # 反向排序

fruits.reverse()     # 反轉順序
```

### 5️⃣ 其他實用操作

```python
# 取得長度
len(fruits)

# 判斷是否存在
if "apple" in fruits:
    print("Yes!")

# 複製 list
copy_fruits = fruits.copy()

# 巢狀列表（Nested List）
matrix = [[1, 2], [3, 4], [5, 6]]
print(matrix[1][0])  # 3
```

---

## 🧠 進階技巧

### 🔸 List Comprehension（列表生成式）

用於快速建立新列表：

```python
squares = [x**2 for x in range(5)]
print(squares)  # [0, 1, 4, 9, 16]
```

### 🔸 遍歷列表

```python
for fruit in fruits:
    print(fruit)
```

### 🔸 結合 enumerate()

```python
for i, fruit in enumerate(fruits):
    print(i, fruit)
```

---

## 🧩 小結

| 操作 | 方法                                 | 範例                        |
| -- | ---------------------------------- | ------------------------- |
| 新增 | `append()`, `insert()`, `extend()` | `fruits.append('kiwi')`   |
| 刪除 | `remove()`, `pop()`, `clear()`     | `fruits.pop(1)`           |
| 搜尋 | `index()`, `count()`               | `nums.index(2)`           |
| 排序 | `sort()`, `reverse()`              | `nums.sort()`             |
| 複製 | `copy()`                           | `copy_list = nums.copy()` |

---

📚 **結論：**
Python 的 `list` 是彈性極高的資料結構，能有效進行資料儲存、操作與轉換。熟練掌握各種方法將顯著提升資料處理與演算法效率。
