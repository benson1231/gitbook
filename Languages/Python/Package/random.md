# 🧩 Python `random` 模組教學

`random` 模組用於在 Python 中生成隨機數、打亂資料或進行隨機選擇。它常用於模擬、遊戲、資料抽樣與機器學習的隨機初始化等情境。

---

## 🔹 1. 匯入模組

```python
import random
```

---

## 🔹 2. 生成隨機數

| 函式                                    | 說明                    | 範例                               |
| ------------------------------------- | --------------------- | -------------------------------- |
| `random.random()`                     | 回傳 0~1 之間的浮點數         | `random.random()` → 0.736        |
| `random.uniform(a, b)`                | 回傳 a~b 之間的浮點數         | `random.uniform(10, 20)` → 14.38 |
| `random.randint(a, b)`                | 回傳 a~b 之間的整數（包含 a, b） | `random.randint(1, 6)` → 3       |
| `random.randrange(start, stop, step)` | 類似 `range()`，隨機取一個數   | `random.randrange(0, 10, 2)` → 6 |

---

## 🔹 3. 從序列中抽樣

| 函式                         | 說明               | 範例                                                 |
| -------------------------- | ---------------- | -------------------------------------------------- |
| `random.choice(seq)`       | 從序列中隨機取一個元素      | `random.choice(['A', 'B', 'C'])` → 'B'             |
| `random.choices(seq, k=n)` | 從序列中取 n 個元素（可重複） | `random.choices([1, 2, 3], k=5)` → [2, 1, 3, 2, 1] |
| `random.sample(seq, k=n)`  | 從序列中取 n 個元素（不重複） | `random.sample(range(10), 3)` → [7, 1, 4]          |
| `random.shuffle(seq)`      | 原地打亂序列順序         | `lst = [1,2,3]; random.shuffle(lst)`               |

---

## 🔹 4. 控制隨機性（設定種子）

使用 `random.seed()` 可讓隨機結果可重現。

```python
random.seed(42)
print(random.random())  # 每次都會輸出相同結果
```

應用場景：

* 機器學習模型初始化
* 測試時需要固定隨機結果

---

## 🔹 5. 常見應用範例

### 🎲 模擬擲骰子

```python
import random
for _ in range(5):
    print(random.randint(1, 6))
```

### 🧠 隨機密碼生成

```python
import random, string
password = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
print(password)
```

### 📊 隨機抽樣分析

```python
data = list(range(1, 101))
sample = random.sample(data, 10)
print(sample)
```

---

## 🔹 6. 延伸模組：`secrets` 與 `numpy.random`

* **`secrets`**：用於產生更安全的隨機數（密碼、token）。
* **`numpy.random`**：支援多種分佈（常態分佈、均勻分佈等），常用於科學運算與機器學習。

---

## ✅ 總結

| 功能      | 函式                                               |
| ------- | ------------------------------------------------ |
| 生成隨機浮點數 | `random.random()`, `random.uniform()`            |
| 生成隨機整數  | `random.randint()`, `random.randrange()`         |
| 抽樣與打亂   | `choice()`, `choices()`, `sample()`, `shuffle()` |
| 固定隨機結果  | `random.seed()`                                  |

---

📘 **延伸閱讀**：

* [Python 官方文件：random 模組](https://docs.python.org/3/library/random.html)
* [Python secrets 模組](https://docs.python.org/3/library/secrets.html)
* [NumPy random 子模組](https://numpy.org/doc/stable/reference/random/index.html)
