# Higher-order Functions and Lambda Function

在 Python 中，有幾個常用的「高階函數（Higher-order Functions）」可用於以簡潔方式處理資料集合。這些函數通常搭配 Lambda（匿名函數）使用，能夠提升程式的表達力與可讀性。

---

## 一、Lambda 函式（匿名函數）

Lambda 是一種不需命名即可定義的簡單函式，常用於一次性處理資料。

### 語法：

```python
lambda 參數: 運算式
```

### 範例：

```python
square = lambda x: x ** 2
print(square(5))  # 輸出 25
```

---

## 二、map()

`map()` 對序列中每個元素套用指定函式，並回傳新的可迭代物件。

### 語法：

```python
map(function, iterable)
```

### 範例：

```python
nums = [1, 2, 3, 4]
squares = list(map(lambda x: x ** 2, nums))
print(squares)  # [1, 4, 9, 16]
```

---

## 三、filter()

`filter()` 用來篩選符合條件的元素。

### 語法：

```python
filter(function, iterable)
```

### 範例：

```python
nums = [1, 2, 3, 4, 5]
evens = list(filter(lambda x: x % 2 == 0, nums))
print(evens)  # [2, 4]
```

---

## 四、reduce()

`reduce()` 將序列中的元素依序合併為單一值，需搭配 `functools` 使用。

### 語法：

```python
from functools import reduce
reduce(function, iterable)
```

### 範例：

```python
from functools import reduce
nums = [1, 2, 3, 4]
product = reduce(lambda x, y: x * y, nums)
print(product)  # 24
```

---

## 五、常見用途比較

| 函數     | 功能說明         | 適合情境                    |
| ------ | ------------ | ----------------------- |
| lambda | 建立簡單函式（不具名）  | 適合短小運算如 `lambda x: x+1` |
| map    | 元素逐一轉換       | 對列表資料進行函式套用             |
| filter | 條件篩選         | 擷取符合條件的項目               |
| reduce | 累積處理（如總和、乘積） | 對所有資料進行彙總（如加總、乘積）       |

---

這些高階函數可與 Python 的列表推導式或資料分析工具（如 Pandas）搭配使用，在保留程式簡潔性的同時，也能強化程式邏輯的表達力。
