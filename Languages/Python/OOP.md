# Object-Oriented Programming(OOP)

Python 是一種支援物件導向程式設計（OOP, Object-Oriented Programming）的語言，允許開發者透過「類別（class）」來封裝資料與行為，構建模組化、可擴展、可重用的程式結構。以下將深入介紹 OOP 的核心概念與 Python 中的實作方式。

---

## 一、物件導向核心概念（OOP）

| 概念                | 說明                                   |
| ----------------- | ------------------------------------ |
| **Class**         | 類別是一個藍圖，定義了物件的屬性與方法，例如人類、動物、車子等概念抽象。 |
| **Object**        | 類別的實體（Instance），代表一個具體的個體，例如一個特定的人。  |
| **Attribute**     | 物件的狀態，例如姓名、年齡等（變數）。                  |
| **Method**        | 操作屬性的函式（功能），例如計算、顯示資訊等。              |
| **Encapsulation** | 封裝：將屬性與方法綁在一起，限制外部直接存取以保護資料。         |
| **Inheritance**   | 繼承：子類別自動取得父類別屬性與方法，可進一步擴充或修改。        |
| **Polymorphism**  | 多型：相同方法名稱可根據物件不同而展現不同行為。             |

---

## 二、基本語法結構與說明

```python
class Person:
    def __init__(self, name, age):
        self.name = name      # 屬性
        self.age = age

    def greet(self):          # 方法
        print(f"Hello, my name is {self.name} and I am {self.age} years old.")

# 建立物件（實例化）
p = Person("Alice", 30)
p.greet()
```

* `class` 宣告類別名稱，命名通常採用大駝峰命名（如：MyClass）。
* `__init__` 為初始化方法，當呼叫 `Person()` 時自動執行。
* `self` 代表物件本身，是屬性與方法的參照依據。

---

## 三、屬性與方法的操作

```python
print(p.name)      # 屬性讀取
p.age = 31         # 屬性修改
p.greet()          # 方法呼叫
```

### 類別層級屬性

```python
class Circle:
    pi = 3.14159   # 類別屬性，共享所有實例

    def __init__(self, radius):
        self.radius = radius

    def area(self):
        return Circle.pi * self.radius ** 2
```

---

## 四、繼承與方法覆寫（Override）

```python
class Student(Person):
    def __init__(self, name, age, student_id):
        super().__init__(name, age)       # 呼叫父類別建構子
        self.student_id = student_id

    def greet(self):                      # 覆寫父類方法
        print(f"I am student #{self.student_id}, my name is {self.name}.")

s = Student("Bob", 20, "S123")
s.greet()
```

---

## 五、私有屬性與封裝

* `_variable`：慣例上的「保護屬性」，提醒僅供內部使用。
* `__variable`：觸發 Python 名稱改寫（Name Mangling），成為私有屬性。

```python
class Bank:
    def __init__(self):
        self._balance = 0          # protected
        self.__pin = 1234          # private

    def show(self):
        print(f"Balance: {self._balance}, PIN: {self.__pin}")
```

---

## 六、類別與靜態方法（@classmethod, @staticmethod）

```python
class Math:
    @staticmethod
    def add(x, y):
        return x + y

    @classmethod
    def description(cls):
        return f"This is class: {cls.__name__}"
```

* `@staticmethod`：不需要 `self`，用來定義與類別無關的工具函式。
* `@classmethod`：接收 `cls` 參數，適合用於建構替代建構子或讀取類別屬性。

---

Python 的物件導向設計簡潔、強大，是構建複雜應用、模組與架構的基石。熟練 class 的定義與操作，有助於編寫具可讀性與可維護性的程式碼。
