# Python 字串（String）教學

## 📘 基本介紹

在 Python 中，**字串（string）** 是以單引號 `' '` 或雙引號 `" "` 括起的文字資料。

```python
name = 'Alice'
greeting = "Hello, world!"
```

字串是**不可變物件（immutable）**，這表示一旦建立就不能被修改。

---

## 🧩 建立字串

```python
# 單引號或雙引號
text1 = 'Hello'
text2 = "World"

# 三引號可建立多行字串
description = '''This is a
multi-line
string.'''
```

---

## 🔍 索引與切片（Indexing & Slicing）

字串可以像 list 一樣用索引取值。

```python
word = 'Python'
print(word[0])   # P
print(word[-1])  # n
print(word[0:3]) # Pyt
```

> `[:]` 可用於切片，支援步長（step）設定，例如 `word[::2]` 會輸出 `Pto`。

---

## 🔤 常用方法

| 方法                            | 說明        | 範例                                         |
| :---------------------------- | :-------- | :----------------------------------------- |
| `len()`                       | 回傳字串長度    | `len('apple') → 5`                         |
| `upper()`                     | 轉大寫       | `'hello'.upper() → 'HELLO'`                |
| `lower()`                     | 轉小寫       | `'HELLO'.lower() → 'hello'`                |
| `capitalize()`                | 首字母大寫     | `'python'.capitalize() → 'Python'`         |
| `title()`                     | 每個單字首字母大寫 | `'hello world'.title() → 'Hello World'`    |
| `strip()`                     | 移除前後空白    | `' hello '.strip() → 'hello'`              |
| `replace(a,b)`                | 取代子字串     | `'hi tom'.replace('tom','sam') → 'hi sam'` |
| `split()`                     | 拆成 list   | `'a,b,c'.split(',') → ['a','b','c']`       |
| `'sep'.join(list)`            | 串接字串      | `'-'.join(['a','b']) → 'a-b'`              |
| `startswith()` / `endswith()` | 判斷開頭或結尾   | `'apple'.startswith('a') → True`           |
| `find()` / `index()`          | 尋找子字串位置   | `'banana'.find('na') → 2`                  |

---

## 🧮 格式化字串（String Formatting）

```python
# f-string（推薦）
name = 'Alice'
age = 25
print(f"My name is {name}, and I am {age} years old.")

# format() 方法
print("My name is {}, and I am {}.".format(name, age))

# 百分比格式化
print("My name is %s, and I am %d." % (name, age))
```

---

## 🧠 判斷與檢查

```python
text = 'Python3'

print(text.isalpha())  # False（因含數字）
print('Hello'.isalpha())  # True
print('123'.isdigit())    # True
print('python'.islower()) # True
```

---

## 🪄 字串反轉

```python
s = 'hello'
print(s[::-1])  # 'olleh'
```

---

## 🧰 多行處理與跳脫字元

```python
# 換行與跳脫字元
print('Line1\nLine2')
print('He said: \"Hi!\"')

# 原始字串（不解析跳脫字元）
path = r'C:\\Users\\Alice'
```

---

## 💡 小結

* 字串是不可變的。
* 常用方法能快速處理大小寫、取代與分割。
* `f-string` 是最方便的格式化方式。
* 善用切片與方法組合，可大幅提升文字處理效率。
