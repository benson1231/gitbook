# Dunder

DUNDER（Double UNDERscore）方法是 Python 類別中的特殊方法，格式為 `__method__`。這些方法讓開發者可以定義類別在特定操作下的行為，例如建立物件、加總、比較、列印等，是物件導向設計的核心。

---

## 一、`__init__`：初始化物件屬性

```python
class Circle:
    def __init__(self, radius):
        self.radius = radius
```

* 建立物件時會自動執行 `__init__`，用於設定初始值。
* `self` 是代表目前的物件實例。

---

## 二、`__str__`：物件的可讀字串表示

```python
    def __str__(self):
        return f"Circle with radius {self.radius}"
```

* 透過 `print(obj)` 或 `str(obj)` 會呼叫此方法。
* 用來提供對使用者友善的輸出格式。

---

## 三、`__repr__`：物件的正式表示（供開發者與除錯使用）

```python
    def __repr__(self):
        return f"Circle({self.radius})"
```

* 透過 `repr(obj)` 或互動式環境呼叫時觸發。
* 通常提供能夠重建物件的語法。

---

## 四、其他常見 DUNDER 方法

| 方法名稱          | 用途說明                  | 範例操作          |
| ------------- | --------------------- | ------------- |
| `__len__`     | 定義 `len(obj)` 行為      | `len(my_obj)` |
| `__eq__`      | 定義 `==` 等值比較邏輯        | `a == b`      |
| `__lt__`      | 定義 `<` 小於比較邏輯         | `a < b`       |
| `__add__`     | 定義 `+` 加法運算邏輯         | `a + b`       |
| `__getitem__` | 支援索引操作                | `obj[0]`      |
| `__call__`    | 讓物件本身可像函式一樣呼叫         | `obj()`       |
| `__del__`     | 當物件被刪除時觸發，用於資源釋放（不常用） | `del obj`     |

---

## 五、完整類別範例

```python
class Circle:
    def __init__(self, radius):
        self.radius = radius

    def __str__(self):
        return f"圓的半徑為 {self.radius}"

    def __repr__(self):
        return f"Circle({self.radius})"

    def __eq__(self, other):
        return self.radius == other.radius

    def __add__(self, other):
        return Circle(self.radius + other.radius)
```

```python
c1 = Circle(3)
c2 = Circle(4)
print(c1 + c2)  # 圓的半徑為 7
print(c1 == c2) # False
```

---

DUNDER 方法是 Python 強大且彈性的設計元素，能讓自定義類別的行為與內建型態一致，便於開發、除錯與重用。
