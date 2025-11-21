# Distance Metrics

在機器學習與資料分析中，距離衡量方法（Distance Metrics）用於計算樣本間的相似程度。不同的距離定義會對模型結果產生重大影響，尤其在聚類、分類、推薦系統中常被使用。

---

## 1. 歐氏距離（Euclidean Distance）

最常見的距離計算方式，即「直線距離」，適用於連續型數值資料。

### 定義：

$$
d(x, y) = \sqrt{\sum_{i=1}^{n}(x_i - y_i)^2}
$$

### Python 範例：

```python
from scipy.spatial.distance import euclidean
x = [1, 2]
y = [4, 6]
print(euclidean(x, y))  # 輸出: 5.0
```

---

## 2. 曼哈頓距離（Manhattan Distance）

又稱「城市街區距離」，計算各維度差的總和。

### 定義：

$$
d(x, y) = \sum_{i=1}^{n} |x_i - y_i|
$$

### Python 範例：

```python
from scipy.spatial.distance import cityblock
print(cityblock(x, y))  # 輸出: 7
```

---

## 3. 餘弦相似度／距離（Cosine Similarity / Distance）

常用於文字向量或高維稀疏資料，衡量向量方向差異而非長度。

### 相似度定義：

$$
\cos(\theta) = \frac{x \cdot y}{\|x\| \cdot \|y\|}
$$

### 餘弦距離：

$$
d(x, y) = 1 - \cos(\theta)
$$

### Python 範例：

```python
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
cos_sim = cosine_similarity([x], [y])
print(cos_sim)  # 越接近 1 表示越相似
```

---

## 4. 切比雪夫距離（Chebyshev Distance）

計算所有維度中最大絕對差值。

### 定義：

$$
d(x, y) = \max_{i} |x_i - y_i|
$$

### Python 範例：

```python
from scipy.spatial.distance import chebyshev
print(chebyshev(x, y))  # 輸出: 4
```

---

## 5. 漢明距離（Hamming Distance）

用於衡量兩個相同長度的字串或二進位向量間有多少個位置不同。

### Python 範例：

```python
from scipy.spatial.distance import hamming
x = [1, 0, 1, 1]
y = [1, 1, 0, 1]
print(hamming(x, y))  # 輸出: 0.5
```

---

## 6. 距離選擇建議

| 應用場景    | 建議距離     |
| ------- | -------- |
| 數值型特徵   | 歐氏、曼哈頓距離 |
| 類別型特徵   | 漢明距離     |
| 文字／嵌入向量 | 餘弦距離     |
| 噪音敏感任務  | 曼哈頓、切比雪夫 |

---

選擇正確的距離衡量方式對於提升模型效能與準確性至關重要。實際應用中可透過交叉驗證測試哪一種距離對目標任務最為合適。
