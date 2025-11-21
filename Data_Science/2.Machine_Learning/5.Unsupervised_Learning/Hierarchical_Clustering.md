# Hierarchical Clustering

Hierarchical Clustering 是一種基於資料間相似度遞迴地建立階層結構的非監督式分群演算法。它不需要預先指定群數，並可用樹狀圖（Dendrogram）呈現群集關係。

---

## 1. 方法分類

### A. Agglomerative（凝聚式，自下而上）
- 初始每個樣本為一個群
- 依照最近距離兩兩合併群
- 一直合併直到成為一大群或達到指定條件

### B. Divisive（分裂式，自上而下）
- 初始為所有樣本一群
- 持續將群拆分直到每個樣本為獨立群
- 計算量高，較少使用

---

## 2. 距離計算方式
| 名稱         | 說明                                |
|--------------|-------------------------------------|
| Euclidean    | 歐式距離（最常用）                 |
| Manhattan    | 絕對距離                           |
| Cosine       | 用於文字或稀疏資料相似度評估       |

## 3. Linkage（群間距離定義）
| 方法         | 群與群之間的距離如何定義            |
|--------------|--------------------------------------|
| Single       | 兩群之間最近距離                    |
| Complete     | 兩群之間最遠距離                    |
| Average      | 所有成對距離的平均值                |
| Ward         | 最小化群內平方誤差總和（推薦）      |

---

## 4. 實作範例（使用 Scikit-learn + scipy）
```python
from sklearn.datasets import load_iris
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

# 載入資料
X, _ = load_iris(return_X_y=True)

# 計算 linkage 矩陣
Z = linkage(X, method='ward')

# 畫出樹狀圖
plt.figure(figsize=(10, 5))
dendrogram(Z)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Sample Index")
plt.ylabel("Distance")
plt.show()
```

---

## 5. 優缺點
### 優點：
- 不需事先指定群數（可視 dendrogram 決定）
- 可產生具階層關係的群分類圖
- 適合中小型資料視覺化

### 缺點：
- 計算與記憶體成本高（O(n²)）
- 對離群點與距離度量敏感
- 不適合大規模資料集

---

## 6. 決定群數的方法
- 觀察 dendrogram 中的垂直切割距離
- 使用 `scipy.cluster.hierarchy.fcluster()` 設定距離或群數閾值分群

---

Hierarchical Clustering 適合探索群體結構與資料間的關係，特別在資料量不大、需要階層式解釋時效果極佳。
