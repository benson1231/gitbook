# Feature Encoding

在機器學習中，資料通常需要轉換成數值形式才能進行模型訓練。這一過程稱為「特徵編碼（Feature Encoding）」，其中常見的方式包括 One-Hot Encoding、Label Encoding、Embedding，以及特徵縮放（Feature Scaling）等方法。

---

## 1. One-Hot Encoding（一熱編碼）

將類別變數轉換為稀疏的二進位向量，只有對應類別的位置為 1，其餘皆為 0。

### 範例：

類別：`['cat', 'dog', 'bird']`

```text
cat  → [1, 0, 0]
dog  → [0, 1, 0]
bird → [0, 0, 1]
```

### 優點：

* 無序性（適合沒有內在順序的分類）

### 缺點：

* 維度爆炸：若類別過多會導致維度過高。

### TensorFlow 實作：

```python
import tensorflow as tf
from sklearn.preprocessing import OneHotEncoder
import numpy as np

data = np.array(['cat', 'dog', 'bird']).reshape(-1, 1)
ohe = OneHotEncoder(sparse=False)
encoded = ohe.fit_transform(data)
print(encoded)
```

---

## 2. Label Encoding（標籤編碼）

將每個類別編碼為一個整數值。

### 範例：

類別：`['cat', 'dog', 'bird']`

```text
cat  → 0
dog  → 1
bird → 2
```

### 優點：

* 節省記憶體空間

### 缺點：

* 模型可能會誤認為類別間有順序關係，不適用於無序分類特徵

### Python 實作：

```python
from sklearn.preprocessing import LabelEncoder

data = ['cat', 'dog', 'bird']
le = LabelEncoder()
labels = le.fit_transform(data)
print(labels)
```

---

## 3. Embedding（嵌入表示）

用於將類別或詞語轉換為稠密的向量空間，適用於大量類別（如文字、物件 ID 等）

### 特點：

* 經由訓練學習出最佳的低維表示
* 常用於 NLP 與深度學習模型中

### TensorFlow 實作：

```python
embedding_layer = tf.keras.layers.Embedding(input_dim=1000, output_dim=64)
```

---

## 4. 特徵縮放（Feature Scaling）

特徵縮放將不同範圍的數值標準化，使其落在統一範圍內，以提升模型收斂速度與效能。常見方法如下：

### (1) Min-Max Scaling（最小最大縮放）

將特徵縮放到 \[0, 1] 範圍。

```python
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
scaled_X = scaler.fit_transform(X)
```

### (2) Standardization（標準化）

將特徵轉換為平均值為 0，標準差為 1 的分布。

```python
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
standardized_X = scaler.fit_transform(X)
```

### 選擇依據：

| 方法              | 特徵範圍控制 | 對離群值敏感性 | 適用場景             |
| --------------- | ------ | ------- | ---------------- |
| Min-Max Scaling | 是      | 高       | 深度學習、距離相關模型（KNN） |
| Standardization | 否      | 低       | 線性模型、邏輯回歸、SVM 等  |

---

## 總結

| 編碼方式            | 優點             | 缺點          | 適用情境        |
| --------------- | -------------- | ----------- | ----------- |
| One-Hot         | 簡單直觀、無序        | 維度高、稀疏矩陣    | 少量類別的分類特徵   |
| Label Encoding  | 計算效率高、節省空間     | 可能誤導模型有序性   | 有明確順序的類別    |
| Embedding       | 支援大量類別、能學習語意關係 | 需訓練模型、解釋性較差 | 自然語言、ID 特徵等 |
| Feature Scaling | 提升訓練效率、消除單位差異  | 需搭配模型特性使用   | 所有數值型特徵     |

---

良好的特徵編碼與縮放策略是機器學習成功的關鍵，根據資料特性選擇適當方法有助於模型效能最大化。

