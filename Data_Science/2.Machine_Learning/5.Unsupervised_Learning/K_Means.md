# K-Means

K-Means 是一種非監督式學習演算法，用於將資料分為 K 個群集（clusters）。它依據資料點之間的距離，將資料分配到距離最近的中心點所屬群集，並持續更新中心直到收斂。

---

## 1. 演算法流程

1. 隨機選擇 K 個初始中心點（centroids）。
2. 將每個資料點指派給距離最近的中心點形成群集。
3. 計算每個群集內資料點的平均值作為新的中心點。
4. 重複步驟 2-3，直到群集不再變動或達到最大迭代次數。

---

## 2. 優缺點

| 優點         | 缺點                |
| ---------- | ----------------- |
| 計算速度快、易於實作 | 需事先指定 K 值         |
| 適合大量資料     | 對初始值敏感，可能陷入局部最小值  |
| 群集結果易解釋    | 假設資料為凸型，對非線性邊界不友善 |

---

## 3. Python 實作

```python
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
import matplotlib.pyplot as plt

# 產生模擬資料
X, _ = make_blobs(n_samples=300, centers=4, cluster_std=0.6, random_state=0)

# 建立 KMeans 模型
kmeans = KMeans(n_clusters=4, random_state=0)
kmeans.fit(X)

# 視覺化分群結果
plt.scatter(X[:, 0], X[:, 1], c=kmeans.labels_, cmap='viridis')
plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], s=200, c='red')
plt.title("K-Means Clustering")
plt.show()
```

---

## 4. 如何選擇最佳 K 值

可透過 **肘部法則（Elbow Method）** 選定最佳 K 值：

1. 計算不同 K 值下的總內部平方誤差（Inertia）。
2. 繪製 K 對應誤差圖形。
3. 找到誤差下降趨勢明顯變緩的 "肘部"，該 K 即為最佳選擇。

```python
inertia = []
for k in range(1, 10):
    kmeans = KMeans(n_clusters=k).fit(X)
    inertia.append(kmeans.inertia_)

plt.plot(range(1, 10), inertia, marker='o')
plt.title("Elbow Method")
plt.xlabel("Number of Clusters")
plt.ylabel("Inertia")
plt.show()
```

---

## 5. 應用範例

* 顧客分群（Customer Segmentation）
* 影像壓縮（Color Quantization）
* 文件主題分群（Text Clustering）
* 社群偵測與市場分析

---

K-Means 是最常見的分群方法之一，適用於大多數線性可分群的資料分析情境。雖然對初始值與 K 值選擇敏感，但搭配肘部法與多次初始化可有效提升準確性與穩定性。
