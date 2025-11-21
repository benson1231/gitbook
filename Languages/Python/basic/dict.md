# 🧩 Python 字典（Dictionary）教學

## 📘 基本介紹

`dict`（字典）是 Python 中用於儲存「鍵值對（key-value pairs）」的資料結構。每個鍵（key）都是唯一的，用來對應一個值（value）。字典屬於可變（mutable）型別，因此可以新增、修改或刪除項目。

```python
# 建立字典
person = {
    "name": "Alice",
    "age": 25,
    "city": "Taipei"
}
```

---

## 🔍 常用方法

### 1. `get()` — 取得值（避免錯誤）

```python
print(person.get("name"))       # 輸出: Alice
print(person.get("country", "N/A"))  # 若 key 不存在，回傳預設值 N/A
```

### 2. `keys()`、`values()`、`items()`

```python
print(person.keys())   # dict_keys(['name', 'age', 'city'])
print(person.values()) # dict_values(['Alice', 25, 'Taipei'])
print(person.items())  # dict_items([('name', 'Alice'), ('age', 25), ('city', 'Taipei')])
```

### 3. `pop()` — 移除指定鍵並回傳其值

```python
age = person.pop("age")
print(age)     # 輸出: 25
print(person)  # {'name': 'Alice', 'city': 'Taipei'}
```

### 4. `popitem()` — 移除最後一組鍵值對

```python
last = person.popitem()
print(last)    # 例如 ('city', 'Taipei')
```

### 5. `update()` — 合併或更新字典內容

```python
person.update({"city": "Tokyo", "job": "Engineer"})
print(person)
# 輸出: {'name': 'Alice', 'city': 'Tokyo', 'job': 'Engineer'}
```

### 6. `clear()` — 清空字典

```python
person.clear()
print(person)  # 輸出: {}
```

---

## 🔁 迴圈遍歷字典

```python
for key, value in person.items():
    print(f"{key}: {value}")
```

---

## ⚙️ 字典生成式（Dictionary Comprehension）

```python
squares = {x: x**2 for x in range(5)}
print(squares)  # {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}
```

## 🔸 將兩個列表合併成字典

```python
keys = ["name", "age", "city"]
values = ["Alice", 25, "Taipei"]
user = dict(zip(keys, values))
```

---

## 🧠 小結

| 操作    | 方法                                | 說明             |
| ----- | --------------------------------- | -------------- |
| 取得值   | `get(key, default)`               | 取值且避免 KeyError |
| 新增/更新 | `update()`                        | 合併或更新內容        |
| 刪除    | `pop(key)` / `popitem()`          | 移除項目           |
| 清空    | `clear()`                         | 清除所有項目         |
| 取鍵值   | `keys()` / `values()` / `items()` | 取得所有鍵、值或鍵值對    |
