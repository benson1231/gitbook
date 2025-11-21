# Normalization

在資料前處理中，Normalization（正規化）是常見的技術，用於將不同特徵的數值範圍轉換至相同尺度，避免模型受到特徵數值大小影響。

---

## 1. 為什麼需要正規化？

* 某些演算法（如 KNN、SVM、線性回歸、神經網路）對數值尺度敏感
* 不同單位的特徵會導致模型偏向數值較大的特徵
* 使得梯度下降更穩定、收斂更快

---

## 2. 常見的正規化方法

### (1) Min-Max Normalization（最小最大正規化）

將資料壓縮到固定範圍（通常是 0 到 1）。

$$
x_{norm} = \frac{x - x_{min}}{x_{max} - x_{min}}
$$

```python
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X)
```

### (2) Z-score Standardization（標準化）

將資料轉換為平均值為 0，標準差為 1 的分布。

$$
z = \frac{x - \mu}{\sigma}
$$

```python
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_standardized = scaler.fit_transform(X)
```

### (3) Max Abs Scaling（最大絕對值縮放）

適用於稀疏資料（如文字向量），將值壓縮到 -1 到 1 之間。

```python
from sklearn.preprocessing import MaxAbsScaler
scaler = MaxAbsScaler()
X_maxabs = scaler.fit_transform(X)
```

### (4) Robust Scaling（穩健縮放）

使用中位數與四分位距，對離群值不敏感。

```python
from sklearn.preprocessing import RobustScaler
scaler = RobustScaler()
X_robust = scaler.fit_transform(X)
```

---

## 3. 選擇建議

| 方法               | 適用情境               |
| ---------------- | ------------------ |
| Min-Max Scaling  | 無離群值，特徵需壓縮至特定範圍    |
| Z-score Standard | 常態分布假設適用、多數 ML 演算法 |
| MaxAbs           | 稀疏資料（如文字 TF-IDF）   |
| Robust           | 含離群值資料             |

---

## 4. 注意事項

* 訓練與測試集應使用相同 scaler
* 正規化後不再保有原始物理意義（如金額、溫度）
* 某些演算法（如樹模型）不需要正規化

---

正規化是有效提升機器學習模型穩定性與準確度的前處理步驟，應根據資料分佈與任務選擇合適的方法。
