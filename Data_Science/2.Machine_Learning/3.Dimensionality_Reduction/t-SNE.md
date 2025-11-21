# t-SNE (t-Distributed Stochastic Neighbor Embedding)

t-SNE 是一種常用的非線性降維技術，特別適合用於高維資料的視覺化，例如將資料從數百維壓縮到 2D 或 3D 空間以進行視覺分析。

---

## 1. 主要用途
- 高維資料的視覺化（常見於圖片、基因、NLP）
- 發現資料中的群集或隱含結構
- 搭配 PCA 作為預處理（加速與穩定）

---

## 2. 原理概念（簡化版）
- 在高維空間中，t-SNE 會衡量樣本間的機率相似度（近的點有較高機率被當成鄰居）
- 在低維空間中重新配置點，使其保有高維空間的局部鄰近關係
- 損失函數使用 Kullback-Leibler divergence 來衡量這種相似度的差異
- 使用梯度下降方式進行最佳化

---

## 3. 主要參數
| 參數         | 說明                                                                 |
|--------------|----------------------------------------------------------------------|
| `perplexity` | 控制鄰近點數量（通常設為 5~50），類似於 kNN 的 k 值                 |
| `n_iter`     | 最大迭代次數（通常設為 1000~5000）                                    |
| `learning_rate` | 學習率（建議介於 10 到 1000 間）                                       |
| `n_components` | 降維後的維度（通常為 2 或 3）                                           |

---

## 4. 範例程式碼（使用 scikit-learn）
```python
from sklearn.manifold import TSNE
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt

# 載入資料
X, y = load_iris(return_X_y=True)

# 建立模型
tsne = TSNE(n_components=2, perplexity=30, learning_rate=200, n_iter=1000)
X_embedded = tsne.fit_transform(X)

# 視覺化
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y)
plt.title("t-SNE Visualization of Iris Dataset")
plt.xlabel("Dimension 1")
plt.ylabel("Dimension 2")
plt.show()
```

---

## 5. 優缺點
### 優點：
- 適合複雜資料的視覺化
- 保留局部結構，比 PCA 更能凸顯群集

### 缺點：
- 不適合用於資料壓縮後進行機器學習建模
- 非常耗時，無法應對大量樣本（>10,000）
- 不具可解釋性（結果無法反推）
- 不保留全域結構（如距離比例）

---

## 6. 常見用途場景
- 圖像資料：MNIST 手寫數字視覺化
- 生物資訊：RNA-seq 或 scRNA-seq 資料群集視覺化
- NLP：word embeddings 的分布觀察（如 Word2Vec）
