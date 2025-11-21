# Hyperparameters & Fine-tuning

在機器學習中，**超參數（Hyperparameters）** 是在模型訓練開始前所設定的參數，與資料本身無直接關係，通常無法透過訓練資料自動學習；而 **微調（Fine-tuning）** 是針對模型訓練過程中進行的優化策略，以提升泛化能力或預測效果。

---

### 1. 超參數與參數的區別

| 類型   | 說明                  | 設定方式            |
| ---- | ------------------- | --------------- |
| 模型參數 | 在訓練過程中學習出來的值（如迴歸係數） | 由模型自動學習         |
| 超參數  | 需在訓練前設定，控制學習流程與結構   | 需透過手動、網格搜尋等方式設定 |

---

### 2. 常見的超參數

以下為不同模型中常見的超參數範例：

#### 2.1 線性回歸（含正則化）

* `alpha`：正則化強度（Lasso、Ridge）
* `fit_intercept`：是否包含截距

#### 2.2 決策樹 / 隨機森林

* `max_depth`：樹的最大深度
* `min_samples_split`：節點拆分的最小樣本數
* `n_estimators`：樹的數量（隨機森林、梯度提升樹）

#### 2.3 神經網路

* `learning_rate`：學習率
* `batch_size`：每次更新所使用的樣本數
* `epochs`：訓練週期數
* `optimizer`：優化器選擇（如 SGD, Adam）

---

### 3. 模型微調策略（Fine-tuning Strategies）

微調的目標為找到最佳的超參數組合，使模型能達到最佳效能。

#### 3.1 手動微調

根據經驗與結果觀察進行調整。

#### 3.2 Grid Search（網格搜尋）

列舉所有可能的超參數組合，針對每一組進行交叉驗證。

```python
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import Ridge

param_grid = {'alpha': [0.01, 0.1, 1, 10]}
model = Ridge()
grid = GridSearchCV(model, param_grid, cv=5)
grid.fit(X, y)

print("最佳 alpha：", grid.best_params_)
```

#### 3.3 Randomized Search

隨機從參數空間抽樣若干組合進行評估，效率較高。

#### 3.4 Bayesian Optimization / HyperOpt / Optuna

進階自動搜尋方法，能更有效找到最佳參數。

#### 3.5 微調的三種方式

| 微調方式                          | 說明                             | 適用情境              |
| ----------------------------- | ------------------------------ | ----------------- |
| **凍結預訓練權重（Freeze）**           | 保留預訓練模型的所有權重，只訓練最後幾層（通常是分類層）   | 資料量小，特徵與原始任務相近    |
| **部分微調（Partial Fine-tuning）** | 解凍部分層（通常是靠近輸出端的幾層），同時訓練這些層與分類器 | 資料中等、任務與原始不同但特徵相近 |
| **全模型微調（Full Fine-tuning）**   | 解凍所有層，讓模型在新任務上進行完整再訓練          | 資料量大，或原任務與新任務差異大  |

---

### 4. 過擬合與泛化能力

選擇過大的模型或錯誤的超參數組合，可能導致過擬合（overfitting），即模型在訓練資料上表現良好但泛化能力差。適當微調超參數（如正則化強度、模型深度）是改善過擬合的關鍵策略之一。

---

### 5. 小結

* 超參數需在訓練前設定，會影響模型結構與訓練流程
* 微調是模型訓練中的重要步驟，可提升模型效能與穩定性
* 善用 Grid Search、交叉驗證與正則化技巧可達到良好效果

接下來可考慮學習 AutoML 與 Hyperparameter Optimization 套件（如 Optuna）以進一步提升效率。
