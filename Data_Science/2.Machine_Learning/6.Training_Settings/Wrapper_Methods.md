# Wrapper Methods

Wrapper 方法是一類特徵選擇技術，透過實際訓練機器學習模型來評估特徵子集的效能。與濾波式方法（Filter）不同，Wrapper 方法將特徵選擇視為模型訓練的一部分。

---

### 1. 方法原理

* 給定一個學習演算法（如決策樹、SVM），
* 透過系統性地嘗試不同的特徵子集（例如透過前向選擇或後向刪除），
* 評估每個子集的模型效能（如交叉驗證準確率），
* 選擇最佳子集。

---

### 2. 常見策略

以下展示每一種策略的簡單 Python 範例：

#### a. 前向選擇（Forward Selection）

* 初始特徵集合為空。
* 每一步中，從尚未被選擇的特徵中，選出一個能使模型效能（如交叉驗證準確率）提升最多的特徵。
* 重複進行直到新增任何特徵皆無明顯效能提升，或達到指定的特徵數量。

**Python 範例：使用 SequentialFeatureSelector（sklearn）**

```python
from sklearn.datasets import load_iris
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.linear_model import LogisticRegression

X, y = load_iris(return_X_y=True)
model = LogisticRegression(max_iter=200)
sfs = SequentialFeatureSelector(model, n_features_to_select=2, direction='forward')
sfs.fit(X, y)
print("選中的特徵索引：", sfs.get_support())
```

- 每一步中，從尚未被選擇的特徵中，選出一個能使模型效能（如交叉驗證準確率）提升最多的特徵。
- 重複進行直到新增任何特徵皆無明顯效能提升，或達到指定的特徵數量。

#### b. 後向刪除（Backward Elimination）
- 初始特徵集合包含所有特徵。
- 每一步中，移除一個對模型貢獻最小的特徵（如 p 值最大、重要性最小者）。
- 持續移除直到效能顯著下降，或特徵數目達到設定門檻。

**Python 範例：使用 statsmodels 搭配 p 值**

```python
import statsmodels.api as sm
import pandas as pd
from sklearn.datasets import load_boston

X, y = load_boston(return_X_y=True, as_frame=True)
X = sm.add_constant(X)
model = sm.OLS(y, X).fit()
print(model.summary())
# 根據 p 值手動刪除不顯著特徵
```

- 每一步中，移除一個對模型貢獻最小的特徵（如 p 值最大、重要性最小者）。
- 持續移除直到效能顯著下降，或特徵數目達到設定門檻。

#### c. 遞迴特徵消除（Recursive Feature Elimination, RFE）
- 利用模型（如線性回歸、SVM）給定的特徵重要性，反覆移除最不重要的特徵。
- 每一輪重新訓練模型與更新重要性排名，直到只剩下所需的特徵數。

#### d. 遺傳演算法（Genetic Algorithm）
- 將特徵子集視為基因染色體。
- 利用交配、突變等機制產生新的特徵組合。
- 選擇在每代中表現最好的特徵子集，重複迭代直到收斂或達到指定輪數。

**Python 範例（使用 sklearn-genetic-opt）**
```python
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from geneticalgorithm import geneticalgorithm as ga

# 需將特徵子集視為二進位編碼的解
# 這裡僅展示簡略邏輯，可搭配自定義 fitness 函數
```
- 利用交配、突變等機制產生新的特徵組合。
- 選擇在每代中表現最好的特徵子集，重複迭代直到收斂或達到指定輪數。

---

### 3. 優缺點比較

| 優點                                | 缺點                                |
|-------------------------------------|-------------------------------------|
| 與指定模型緊密整合，評估準確       | 計算成本高                          |
| 可發現與模型高度相關的最佳子集     | 容易過擬合（需交叉驗證）            |

---

### 4. Python 範例：遞迴特徵消除（RFE）

```python
from sklearn.datasets import load_iris
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression

X, y = load_iris(return_X_y=True)
model = LogisticRegression(max_iter=200)
rfe = RFE(model, n_features_to_select=2)
rfe.fit(X, y)

print("被選中的特徵索引：", rfe.support_)
print("特徵排名：", rfe.ranking_)
```

---

### 5. 應用情境

* 模型為中心的特徵挑選
* 搭配資料量適中、特徵不多的場景
* 適合用於進行模型壓縮、提升可解釋性

---

Wrapper 方法在特徵子集評估上表現強大，但計算資源消耗高。實務中常搭配 Filter 與 Embedded 方法（如 Lasso）做綜合特徵選擇策略。
