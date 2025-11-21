# Python 物件導向程式設計（OOP）教學

Python 是一個支援物件導向（Object-Oriented Programming, OOP）的語言，允許你以「類別（class）」和「物件（object）」的方式組織程式碼，使程式更具可重用性與維護性。

---

## 一、OOP 基本概念

| 概念                    | 說明                                       |
| --------------------- | ---------------------------------------- |
| **Class（類別）**         | 定義物件的藍圖或模板，包含屬性（attributes）與方法（methods）。 |
| **Object（物件）**        | 由類別實例化（instantiate）而來的個體。                |
| **Attribute（屬性）**     | 用來儲存物件的資料（變數）。                           |
| **Method（方法）**        | 類別中的函式，用來定義物件行為。                         |
| **Inheritance（繼承）**   | 子類別可以繼承父類別的屬性與方法。                        |
| **Encapsulation（封裝）** | 隱藏內部資料與實作細節，只提供必要介面。                     |
| **Polymorphism（多型）**  | 不同類別可以定義相同名稱的方法，但表現不同行為。                 |

---

## 二、定義類別與建立物件

```python
class Dog:
    # 初始化方法（constructor）
    def __init__(self, name, age):
        self.name = name  # 屬性
        self.age = age

    # 實例方法
    def bark(self):
        return f"{self.name} is barking!"

# 建立物件
dog1 = Dog("Buddy", 3)
print(dog1.bark())  # Buddy is barking!
```

### 🔹 `init` 方法

* 是「建構子（constructor）」，在物件建立時自動呼叫。
* 用於初始化物件的屬性。

---

## 三、常見特殊方法（Magic Methods / Dunder Methods）

| 方法           | 說明                 | 範例                         |
| ------------ | ------------------ | -------------------------- |
| `__init__()` | 初始化物件屬性            | `def __init__(self, name)` |
| `__str__()`  | 回傳物件的可讀字串          | `print(obj)` 時呼叫           |
| `__repr__()` | 回傳物件的正式表示（通常供除錯使用） | `repr(obj)`                |
| `__len__()`  | 讓物件可被 `len()` 呼叫   | `len(obj)`                 |
| `__add__()`  | 定義加法運算行為           | `obj1 + obj2`              |
| `__eq__()`   | 定義等號 `==` 的比較方式    | `obj1 == obj2`             |
| `__call__()` | 讓物件可被當作函式呼叫        | `obj()`                    |

### 範例： `repr` 與 `str`

```python
class Dog:
    def __init__(self, name, age):
        self.name = name
        self.age = age

    def __repr__(self):
        return f"Dog(name={self.name!r}, age={self.age!r})"

    def __str__(self):
        return f"{self.name}, {self.age} years old"

dog = Dog("Buddy", 3)
print(str(dog))   # Buddy, 3 years old
print(repr(dog))  # Dog(name='Buddy', age=3)
```

---

## 四、繼承（Inheritance）

```python
class Animal:
    def speak(self):
        return "Some sound"

class Dog(Animal):
    def speak(self):  # 覆寫（Override）父類方法
        return "Woof!"

class Cat(Animal):
    def speak(self):
        return "Meow!"

dog = Dog()
cat = Cat()
print(dog.speak())  # Woof!
print(cat.speak())  # Meow!
```

---

## 五、封裝（Encapsulation）

Python 沒有嚴格的私有屬性概念，但以命名慣例達成：

```python
class Account:
    def __init__(self, owner, balance):
        self.owner = owner
        self.__balance = balance  # 私有屬性（建議外部不要直接存取）

    def deposit(self, amount):
        self.__balance += amount

    def get_balance(self):
        return self.__balance

acct = Account("Alice", 1000)
acct.deposit(500)
print(acct.get_balance())  # 1500
```

---

## 六、多型（Polymorphism）

多型允許不同物件在呼叫相同方法名稱時表現出不同的行為：

```python
for animal in [Dog(), Cat()]:
    print(animal.speak())
```

---

## 七、實務建議

* 使用 `__repr__()` 方便除錯與日誌輸出。
* 使用 `@property` 建立安全的屬性存取介面。
* 將共用邏輯放入父類別，提高程式可維護性。

---

✅ **小結**
Python 的 OOP 架構使程式結構更清晰，能夠有效封裝資料、擴展功能並提高可重用性。熟悉如 `__init__`、`__repr__`、`__str__` 等特殊方法，是撰寫高品質 Python 類別的基礎。
