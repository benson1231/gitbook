# NumPy

`NumPy`（Numerical Python）是 Python 中進行數值運算與陣列操作的核心套件，廣泛應用於資料科學、機器學習與影像處理等領域。

---

## 一、安裝

```bash
pip install numpy
```

---

## 二、主要功能與特性

| 功能分類             | 說明                        |
| ---------------- | ------------------------- |
| 多維陣列（ndarray）    | 高效儲存與操作數值資料               |
| 廣播（broadcasting） | 支援不同形狀陣列之間的運算             |
| 向量化計算            | 加速迴圈運算，支援矩陣加減乘除等          |
| 隨機數生成            | 亂數分布、重複抽樣等                |
| 線性代數運算           | dot、inv、eig、svd 等矩陣運算     |
| 統計函數             | mean、std、percentile 等統計方法 |
| 資料 IO            | 支援文字、CSV、二進位檔案格式的讀寫       |

---

## 三、ndarray 建立與操作

```python
import numpy as np

arr = np.array([[1, 2, 3], [4, 5, 6]])
print(arr.shape)     # (2, 3)
print(arr.dtype)     # int64
```

### 常見建立方式：

```python
np.zeros((3, 3))      # 全零矩陣
np.ones((2, 4))       # 全一矩陣
np.eye(3)             # 單位矩陣
np.arange(0, 10, 2)   # 等差數列
np.linspace(0, 1, 5)  # 等距分割數列
```

---

## 四、索引、切片與邏輯運算

```python
arr[0, 1]           # 單一元素
arr[:, 1]           # 取第 2 欄
arr[arr > 3]        # 條件過濾
```

---

## 五、基本運算與廣播

```python
x = np.array([1, 2, 3])
y = np.array([10])
print(x + y)        # 自動廣播為 [11, 12, 13]
```

---

## 六、矩陣與線性代數

```python
A = np.array([[1, 2], [3, 4]])
B = np.linalg.inv(A)
C = np.dot(A, B)
```

---

## 七、統計與隨機數

```python
np.mean(arr)
np.std(arr)
np.random.randint(0, 10, (3, 3))
np.random.normal(loc=0, scale=1, size=1000)
```

---

## 八、檔案操作

```python
np.savetxt("data.csv", arr, delimiter=",")
data = np.loadtxt("data.csv", delimiter=",")
```

---

## 九、應用場景

* 科學計算（數學模擬、物理建模）
* 機器學習（特徵矩陣處理、模型運算）
* 影像處理（搭配 OpenCV、skimage）
* 財務分析、統計數據運算

---

NumPy 提供高效能、多功能的陣列與矩陣處理介面，是所有數值型 Python 應用的基礎。熟練 NumPy 是進入資料科學與 AI 領域的第一步。
