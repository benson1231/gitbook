# 🧮 Python `decimal` 模組教學

Python 的 `decimal` 模組提供**高精度的小數運算**功能，用於避免浮點數計算中的誤差問題，特別適合用在**金融、科學計算與精度敏感的應用**中。

---

## 📘 一、基本概念

在 Python 中，內建的 `float` 使用二進位浮點表示方式，容易產生微小的誤差：

```python
>>> 0.1 + 0.2
0.30000000000000004
```

使用 `decimal.Decimal` 可以避免這種問題：

```python
from decimal import Decimal

>>> Decimal('0.1') + Decimal('0.2')
Decimal('0.3')
```

> ⚠️ **注意：** 要以字串建立 Decimal，而不是 float，否則仍會帶入浮點誤差。

---

## ⚙️ 二、建立 Decimal 物件

```python
from decimal import Decimal

# 從字串建立
x = Decimal('10.25')

# 從整數建立
y = Decimal(3)

# 從浮點數建立（不建議）
z = Decimal(0.1)  # ⚠️ 會帶入誤差
```

---

## 🔢 三、基本運算

`Decimal` 支援加減乘除與比較運算，且保留精確值：

```python
from decimal import Decimal

a = Decimal('1.1')
b = Decimal('2.2')

print(a + b)  # 3.3
print(a * b)  # 2.42
print(b / a)  # 2
print(a < b)  # True
```

---

## ⚖️ 四、設定精度與捨入模式

使用 `getcontext()` 調整全域運算設定：

```python
from decimal import Decimal, getcontext

getcontext().prec = 4  # 設定小數精度
getcontext().rounding = 'ROUND_HALF_UP'  # 四捨五入

x = Decimal('1') / Decimal('3')
print(x)  # Decimal('0.3333')
```

常用捨入模式：

* `ROUND_HALF_UP`：四捨五入（常見）
* `ROUND_DOWN`：無條件捨去
* `ROUND_UP`：無條件進位
* `ROUND_HALF_EVEN`：銀行家捨入法

---

## 🧰 五、常用方法

| 方法            | 功能          | 範例                               |
| ------------- | ----------- | -------------------------------- |
| `.quantize()` | 設定小數位數與捨入方式 | `x.quantize(Decimal('0.01'))`    |
| `.sqrt()`     | 平方根         | `Decimal('16').sqrt()` → `4`     |
| `.exp()`      | e 的次方       | `Decimal(1).exp()` → `2.7182...` |
| `.ln()`       | 自然對數        | `Decimal('10').ln()`             |
| `.log10()`    | 以 10 為底的對數  | `Decimal('100').log10()`         |

---

## 💡 六、應用範例：貨幣計算

```python
from decimal import Decimal, getcontext

getcontext().prec = 10

price = Decimal('19.99')
quantity = Decimal('3')
tax_rate = Decimal('0.05')

subtotal = price * quantity

total = subtotal * (Decimal('1') + tax_rate)
print(f"Total: ${total.quantize(Decimal('0.01'))}")  # Total: $62.97
```

---

## 🧩 七、與 float 的比較

| 特性 | `float`     | `Decimal`     |
| -- | ----------- | ------------- |
| 精度 | 有誤差（約 15 位） | 高可控，可設定精度     |
| 速度 | 快           | 稍慢            |
| 用途 | 一般數值運算      | 金融、科學計算、需精度應用 |

---

## ✅ 八、重點總結

* 使用 `Decimal('0.1')` 而非 `Decimal(0.1)` 以避免誤差。
* 可透過 `getcontext()` 設定精度與捨入方式。
* 適合金融、統計、會計與科學精度運算。
* 雖然比 `float` 慢，但提供可預測且一致的結果。

---

📚 **延伸閱讀：**

* 官方文件：[https://docs.python.org/3/library/decimal.html](https://docs.python.org/3/library/decimal.html)
