# UMAP (Uniform Manifold Approximation and Projection)

UMAP 是一種用於非線性降維的演算法，與 t-SNE 類似，但通常更快且更能保留全域結構。它被廣泛應用於資料視覺化與高維資料探索。

---

## 1. 主要用途
- 將高維資料（如 100 維）壓縮至 2D 或 3D 以利視覺化
- 探索潛在群集或結構
- 常用於圖像、聲音、基因、文本資料等領域

---

## 2. 原理簡介
- 基於 Riemann 幾何與拓撲學
- 在高維空間中構建樣本間的模糊拓撲圖（fuzzy simplicial complex）
- 在低維空間中保留拓撲結構，使鄰近關係盡可能一致
- 透過圖優化方式達成降維結果

---

## 3. 常見參數
| 參數              | 說明                                                                 |
|-------------------|----------------------------------------------------------------------|
| `n_neighbors`      | 控制局部性，越小越強調局部結構（類似 perplexity）                   |
| `n_components`     | 降維後的維度（通常為 2 或 3）                                          |
| `min_dist`         | 控制低維空間點與點的最小距離，影響群集緊密度                          |
| `metric`           | 距離計算方式（如 'euclidean'、'cosine'）                              |

---

## 4. 範例程式碼（使用 umap-learn）
```python
import umap
from sklearn.datasets import load_digits
import matplotlib.pyplot as plt

# 載入資料
X, y = load_digits(return_X_y=True)

# 建立模型
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, metric='euclidean')
X_embedded = reducer.fit_transform(X)

# 視覺化
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y, cmap='Spectral', s=5)
plt.title("UMAP projection of Digits dataset")
plt.show()
```

---

## 5. 優缺點
### 優點：
- 執行速度快、可處理大規模資料（t-SNE 不易）
- 保留更多全域結構
- 可用於特徵壓縮後進行機器學習

### 缺點：
- 對 `n_neighbors` 與 `min_dist` 參數敏感
- 在特定資料上會產生不同視覺效果（不如 t-SNE 一致）
- 結果仍較難解釋

---

## 6. 應用場景
- 單細胞 RNA-seq 資料可視化
- 大規模圖像特徵（如 CNN embedding）投影
- 與分類器搭配做降維 + 預測管線（如 UMAP + KNN）
