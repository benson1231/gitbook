# K-Nearest Neighbors (KNN)

KNN（K 最鄰近演算法）是一種非參數的監督式學習方法，可用於分類與回歸任務。

---

### 1. 演算法原理

* 給定一個測試點，KNN 根據訓練集中距離該點最近的 K 個鄰居來進行預測。
* 分類：以多數決投票方式決定類別。
* 回歸：取 K 個鄰居的平均值作為預測。

---

### 2. 距離度量

常見的距離計算方法：

* 歐氏距離（Euclidean Distance）：

$$
d(x, x') = \sqrt{\sum_{i=1}^{n} (x_i - x'_i)^2}
$$

* 曼哈頓距離（Manhattan Distance）
* 明氏距離（Minkowski Distance）

---

### 3. 優缺點

| 優點        | 缺點              |
| --------- | --------------- |
| 簡單直觀，容易實作 | 計算成本高（需儲存所有資料）  |
| 可處理非線性邊界  | 對離群值敏感          |
| 無需模型假設    | 高維資料下效能下降（維度災難） |

---

### 4. Python 實作

```python
from sklearn.neighbors import KNeighborsClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# 載入資料集
iris = load_iris()
X_train, X_test, y_train, y_test = train_test_split(
    iris.data, iris.target, test_size=0.3, random_state=42)

# 建立模型並訓練
knn = KNeighborsClassifier(n_neighbors=3)
knn.fit(X_train, y_train)

# 預測與評估
y_pred = knn.predict(X_test)
print("準確率：", accuracy_score(y_test, y_pred))
```

---

### 5. 重要參數

* `n_neighbors`：設定 K 值（鄰居數量）
* `metric`：距離函數，預設為 'minkowski'
* `weights`：'uniform' 或 'distance'

---

### 6. 應用場景

* 醫學診斷（根據症狀分類疾病）
* 推薦系統（找出相似使用者）
* 圖像分類與物件辨識

---

### 7. 可視化決策邊界（2D 示例）

```python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# 取前兩維進行可視化
X = iris.data[:, :2]
y = iris.target
knn2d = KNeighborsClassifier(n_neighbors=5)
knn2d.fit(X, y)

# 建立網格點進行預測
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.1),
                     np.arange(y_min, y_max, 0.1))
Z = knn2d.predict(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)

plt.figure(figsize=(8, 6))
cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
plt.contourf(xx, yy, Z, cmap=cmap_light, alpha=0.6)
plt.scatter(X[:, 0], X[:, 1], c=y, edgecolor='k')
plt.title("KNN Decision Boundary (K=5)")
plt.xlabel("Feature 1")
plt.ylabel("Feature 2")
plt.show()
```

---

KNN 是一種直觀強大的方法，特別適合樣本數不大、類別邊界複雜的問題。然而隨資料量與維度提升，其效能會快速下降，因此通常搭配特徵選擇或降維方法使用。
