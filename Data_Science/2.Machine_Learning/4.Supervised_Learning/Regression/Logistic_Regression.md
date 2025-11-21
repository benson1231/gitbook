# Logistic Regression

羅吉斯回歸是一種常用於二元分類問題的監督式學習演算法，其目標是預測輸入特徵屬於某一類別的機率。與線性回歸不同的是，它的輸出為介於 0 與 1 之間的機率值。

---

## 1. 模型形式

假設有輸入特徵向量 $\mathbf{x} = [x_1, x_2, \dots, x_n]$，對應的模型為：

$h_\theta(x) = \frac{1}{1 + e^{-z}} \quad \text{其中 } z = \theta_0 + \theta_1 x_1 + \dots + \theta_n x_n$

這個函數被稱為 sigmoid 函數，能將實數輸出轉換為 $[0, 1]$ 的機率值。

---

## 2. 預測與決策邊界

預測類別是根據輸出機率與閾值進行分類，通常閾值設為 0.5：

$$
\hat{y} = \begin{cases}
1 & \text{if } h_\theta(x) \geq 0.5 \\
0 & \text{otherwise}
\end{cases}
$$

---

## 3. 損失函數（Cost Function）

邏輯回歸使用的損失函數是對數損失（log loss）：

$$
J(\theta) = - \frac{1}{m} \sum_{i=1}^m \left[ y^{(i)} \log(h_\theta(x^{(i)})) + (1 - y^{(i)}) \log(1 - h_\theta(x^{(i)})) \right]
$$

其中：

* $m$：樣本數
* $y^{(i)}$：第 $i$ 筆資料的實際類別（0 或 1）
* $h_\theta(x^{(i)})$：模型預測機率

此損失函數對極端預測（如預測為 0 但實際為 1）懲罰較大。

---

## 4. 模型訓練與優化

邏輯回歸常使用梯度下降法（Gradient Descent）來最小化損失函數。

更新公式如下：

$$
\theta_j := \theta_j - \alpha \frac{\partial J(\theta)}{\partial \theta_j}
$$

其中：

* $\alpha$：學習率（learning rate）
* $\frac{\partial J(\theta)}{\partial \theta_j}$：損失函數對參數的偏導數

---

## 5. Python 實作

使用 Scikit-learn 進行邏輯回歸訓練與預測：

```python
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix

# 假設已有特徵 X 與標籤 y
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

model = LogisticRegression()
model.fit(X_train, y_train)

preds = model.predict(X_val)
print("Accuracy:", accuracy_score(y_val, preds))
print("Confusion Matrix:\n", confusion_matrix(y_val, preds))
```

---

## 6. 優點與限制

### 優點：

* 模型簡單，易於實作與解釋
* 適用於二元分類問題
* 可輸出機率解釋結果

### 限制：

* 僅適用於線性可分資料（對複雜非線性資料效果差）
* 容易受離群值影響
* 對於多分類需要擴展為 softmax（多項式邏輯回歸）

---

邏輯回歸雖為簡單模型，但在許多領域（如醫療診斷、信用評分）仍有廣泛應用。
