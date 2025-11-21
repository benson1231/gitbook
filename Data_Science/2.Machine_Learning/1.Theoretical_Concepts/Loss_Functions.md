# Loss Functions

損失函數（Loss Function）是機器學習模型訓練的核心，用來衡量模型預測值與真實值之間的差距。透過最小化損失函數，模型得以不斷調整參數，進而提升預測準確度。

---

## 1. 損失函數 vs 評估指標
| 項目       | 損失函數                           | 評估指標                     |
|------------|------------------------------------|------------------------------|
| 用途       | 模型訓練時最小化的目標             | 測試或驗證模型表現的指標     |
| 常見舉例   | MSE, Cross-Entropy                | Accuracy, F1-score, R²       |
| 是否參與訓練 | 會                                 | 否                           |

---

## 2. 常見損失函數（依任務類型）

### A. 回歸問題

#### 1. Mean Squared Error（MSE）
\[ \text{MSE} = \frac{1}{n} \sum (y_i - \hat{y}_i)^2 \]
- 對大誤差非常敏感（平方懲罰）
- 預設損失於線性回歸、Ridge

#### 2. Mean Absolute Error（MAE）
\[ \text{MAE} = \frac{1}{n} \sum |y_i - \hat{y}_i| \]
- 對離群值較有魯棒性

#### 3. Huber Loss
- 結合 MSE 和 MAE，平衡穩定性與可導性
- 在誤差小時類似 MSE，大時接近 MAE

### B. 分類問題

#### 1. Binary Cross-Entropy（Log Loss）
\[ -[y \log(\hat{y}) + (1 - y) \log(1 - \hat{y})] \]
- 常用於二元分類（如 Logistic Regression）

#### 2. Categorical Cross-Entropy
\[ -\sum y_i \log(\hat{y}_i) \]
- 多類別分類時使用（Softmax 輸出）

#### 3. Hinge Loss
\[ \text{Loss} = \max(0, 1 - y \cdot \hat{y}) \]
- 用於 SVM，鼓勵分類邊界更大

---

## 3. 損失函數選擇原則
| 任務類型   | 損失函數                   |
|------------|----------------------------|
| 回歸       | MSE, MAE, Huber            |
| 二元分類   | Binary Cross-Entropy       |
| 多類別分類 | Categorical Cross-Entropy  |
| 邊界最大化 | Hinge Loss（SVM）         |

---

## 4. 進階延伸（可擴充）
- **Focal Loss**：專門處理類別不平衡的分類問題
- **KL Divergence**：用於機率分布之間的損失衡量（如知識蒸餾）
- **Custom Loss**：深度學習中可自訂損失函數（如重建誤差、對比學習）

---

## 5. 小結
損失函數直接決定模型「學習什麼」，因此選對損失函數等同於正確定義任務目標，是模型訓練成敗的關鍵。
