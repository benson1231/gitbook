# Regularization

Regularization 是一種抑制模型過度擬合（overfitting）的技術，透過在損失函數中加入懲罰項（penalty term），限制模型參數的大小，使模型更加泛化。

---

### 1. 線性模型中的正則化

假設原始的損失函數為均方誤差（MSE）：

$$
\mathcal{L}(\mathbf{w}) = \sum_{i=1}^{n} (y_i - \mathbf{w}^T x_i)^2
$$

加入正則化後的形式：

$$
\mathcal{L}_{\text{reg}}(\mathbf{w}) = \sum_{i=1}^{n} (y_i - \mathbf{w}^T x_i)^2 + \lambda R(\mathbf{w})
$$

其中 \$R(\mathbf{w})\$ 為正則化項，\$\lambda\$ 為控制正則化強度的超參數。

---

### 2. 常見正則化方式

#### L1 正則化（Lasso Regression）

懲罰項為權重的絕對值總和：

$$
R(\mathbf{w}) = \sum_{j=1}^{p} |w_j|
$$

* 具有特徵選擇能力（可將部分係數壓為零）
* 適合特徵數多且稀疏的情境

#### L2 正則化（Ridge Regression）

懲罰項為權重的平方和：

$$
R(\mathbf{w}) = \sum_{j=1}^{p} w_j^2
$$

* 不會產生稀疏解，但能穩定權重
* 適用於多重共線性問題

#### Elastic Net（彈性網）

結合 L1 與 L2：

$$
R(\mathbf{w}) = \alpha \sum |w_j| + (1 - \alpha) \sum w_j^2
$$

* 同時具有特徵選擇與穩定化特性

---

### 3. Python 實作範例

```python
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.datasets import make_regression

# 模擬資料
X, y = make_regression(n_samples=100, n_features=10, noise=10)

# Ridge
ridge = Ridge(alpha=1.0)
ridge.fit(X, y)

# Lasso
lasso = Lasso(alpha=0.1)
lasso.fit(X, y)

# ElasticNet
enet = ElasticNet(alpha=0.1, l1_ratio=0.5)
enet.fit(X, y)
```

---

### 4. 選擇與比較

| 方法          | 懲罰形式    | 優點  | 缺點       |           |           |
| ----------- | ------- | --- | -------- | --------- | --------- |
| L1 (Lasso)  | \$      | w   | \$       | 稀疏性強，特徵選擇 | 在共線性下不穩定  |
| L2 (Ridge)  | \$w^2\$ | 穩定解 | 無法刪除冗餘特徵 |           |           |
| Elastic Net | \$      | w   | + w^2\$  | 綜合兩者優點    | 多一個超參數需調整 |

---

正則化是現代機器學習不可或缺的技術之一，適當地使用能有效提升模型的穩健性與預測效能。
