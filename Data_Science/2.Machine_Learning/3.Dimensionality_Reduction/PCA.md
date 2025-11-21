# PCA

主成分分析（PCA, Principal Component Analysis）是一種常用於降維的線性方法，能將高維資料轉換為較低維的表示，保留大部分變異性。

---

## 一、PCA 原理概述

* **目標**：從原始特徵中找出能解釋資料變異最多的方向（主成分）。
* **特性**：主成分彼此正交（無關）、依變異量排序。
* **用途**：資料視覺化、降維、去除共線性、加速模型訓練。

---

## 二、步驟流程

1. 標準化資料（使平均為 0、變異為 1）
2. 計算共變異數矩陣
3. 求出特徵值與特徵向量
4. 依特徵值排序，選取前 $k$ 個主成分
5. 投影原始資料至新主成分空間

---

## 三、使用 scikit-learn 實作 PCA

### 1. 基本降維與視覺化

```python
from sklearn.decomposition import PCA
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt
import pandas as pd

# 載入資料
iris = load_iris()
X = iris.data
y = iris.target

# 建立 PCA 模型並降為 2 維
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# 繪圖
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=y)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA on Iris Dataset")
plt.show()
```

### 2. 解釋變異比例

```python
print(pca.explained_variance_ratio_)
print("累積貢獻率：", sum(pca.explained_variance_ratio_))
```

---

## 四、非 scikit-learn 傳統實作

以下是純 NumPy 的 PCA 實作，展示其數學背景：

```python
import numpy as np
from sklearn.datasets import load_iris

X = load_iris().data

# 1. 資料標準化
X_meaned = X - np.mean(X, axis=0)

# 2. 計算共變異矩陣
cov_mat = np.cov(X_meaned, rowvar=False)

# 3. 計算特徵值與特徵向量
eigen_vals, eigen_vecs = np.linalg.eigh(cov_mat)

# 4. 排序特徵值與特徵向量（由大至小）
sorted_index = np.argsort(eigen_vals)[::-1]
sorted_eigenvecs = eigen_vecs[:, sorted_index]

# 5. 選取前兩個主成分
n_components = 2
eigenvector_subset = sorted_eigenvecs[:, :n_components]

# 6. 投影至新空間
X_reduced = np.dot(X_meaned, eigenvector_subset)
```

此實作流程與 PCA 原理完全對應，可用於了解數學計算細節或自訂流程。

---

## 五、參數說明

| 參數             | 說明                       |
| -------------- | ------------------------ |
| `n_components` | 降維維度數或保留變異百分比            |
| `whiten=True`  | 將主成分正規化（常用於影像）           |
| `svd_solver`   | 使用的特徵值分解方式（auto, full 等） |

---

## 六、應用場景

* 圖像處理（如臉部辨識、影像壓縮）
* 基因表達資料視覺化與降維
* 資料預處理以解決特徵共線性問題
* 作為其他模型（如 LDA、SVM）的前置降維步驟

---

PCA 是資料科學與機器學習中常見的資料轉換工具之一，理解其數學背景與應用情境有助於更有效地處理高維資料。
