# Training Settings

本章節整理機器學習模型訓練過程中的關鍵設定與技術，這些設定會直接影響模型收斂速度、泛化能力與穩定性。了解這些細節能幫助你提升訓練效率並減少過擬合風險。

---

## 1. 資料切分（Data Splitting）
- **訓練集 / 驗證集 / 測試集**（常見比例如 60/20/20 或 70/15/15）
- 使用 `train_test_split()` 進行隨機切分
- 可加入 `stratify=y` 來維持類別分布一致

```python
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y)
```

---

## 2. 交叉驗證（Cross-Validation）
- 用於小樣本下的模型泛化能力評估
- 常見形式：K-Fold, StratifiedKFold, Leave-One-Out
- `cross_val_score()` 可快速評估模型分數

---

## 3. 批次訓練（Batching）
- 常見於深度學習，將資料分批餵入模型以節省記憶體
- 常見設定：Batch Size = 32, 64, 128
- 與 epoch 數、learning rate 搭配調整效果最佳

---

## 4. 損失函數（Loss Function）
- 衡量模型預測與真實答案的誤差
| 任務類型 | 常見損失函數         |
|----------|----------------------|
| 回歸     | MSE, MAE, Huber      |
| 分類     | Cross-Entropy, Log Loss |

---

## 5. 最佳化演算法（Optimizer）
- 最常見為 Gradient Descent 類方法
- 常見選項：SGD, Adam, RMSProp（深度學習）
- 傳統模型中內部自帶如 L-BFGS, SAG 等演算法

---

## 6. 超參數設定（Hyperparameters）
- 與模型本身設定不同（如：learning rate, tree depth, regularization strength）
- 可透過 Grid Search、Random Search 或 `Bayesian Optimization` 尋找最佳組合

---

## 7. 正則化技巧（Regularization）
- 降低模型過擬合
- L1（Lasso）、L2（Ridge）、Dropout（深度學習）

---

## 8. Early Stopping
- 在驗證集損失不再改善時提早終止訓練，防止過擬合
- 常與模型檢查點（checkpoint）搭配使用

---

## 9. 使用 Pipeline 封裝訓練流程
- 便於特徵處理與模型一體化訓練
- 防止資料洩漏

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('model', LogisticRegression())
])
```

---

Training Settings 的妥善配置能讓模型表現更穩定、實驗更容易重現，是機器學習實務成功的重要基礎。
